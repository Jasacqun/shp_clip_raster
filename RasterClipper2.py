# -*- coding: UTF8 -*-

from osgeo import gdal, gdalnumeric, ogr, gdal_array
from PIL import Image, ImageDraw
import os
import operator
from functools import reduce
import numpy as np
from scipy import misc

gdal.UseExceptions()


# This function will convert the rasterized clipper shapefile
# to a mask for use within GDAL.
def imageToArray(i):
    """
    Converts a Python Imaging Library array to a
    gdalnumeric image.
    """
    a = gdalnumeric.fromstring(i.tobytes(), 'b')
    a.shape = i.im.size[1], i.im.size[0]
    return a


def arrayToImage(a):
    """
    Converts a gdalnumeric array to a
    Python Imaging Library Image.
    """
    i = Image.frombytes('L', (a.shape[1], a.shape[0]), (a.astype('b')).tobytes())
    return i


def world2Pixel(geoMatrix, x, y):
    """
    Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
    the pixel location of a geospatial coordinate
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / xDist)
    return (pixel, line)


#
#  EDIT: this is basically an overloaded
#  version of the gdal_array.OpenArray passing in xoff, yoff explicitly
#  so we can pass these params off to CopyDatasetInfo
#
def OpenArray(array, prototype_ds=None, xoff=0, yoff=0):
    ds = gdal.Open(gdalnumeric.GetArrayFilename(array))

    if ds is not None and prototype_ds is not None:
        if type(prototype_ds).__name__ == 'str':
            prototype_ds = gdal.Open(prototype_ds)
        if prototype_ds is not None:
            gdalnumeric.CopyDatasetInfo(prototype_ds, ds, xoff=xoff, yoff=yoff)
    return ds


def histogram(a, bins=range(0, 256)):
    """
    Histogram function for multi-dimensional array.
    a = array
    bins = range of numbers to match
    """
    fa = a.flat
    n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
    n = gdalnumeric.concatenate([n, [len(fa)]])
    hist = n[1:] - n[:-1]
    return hist

def stretch(a):
  """
  Performs a histogram stretch on a gdalnumeric array image.
  """
  real = np.zeros(a.shape)
  for i in range(len(a)):
      for j in range(len(a[i])):
          val = a[i][j].real
          real[i, j] = val

  hist = histogram(real)
  im = arrayToImage(real)
  lut = []
  for b in range(0, len(hist), 256):
    # step size
    step = reduce(operator.add, hist[b:b+256]) / 255
    # create equalization lookup table
    n = 0
    for i in range(256):
      lut.append(n / step)
      n = n + hist[i+b]
  im = im.point(lut)
  return imageToArray(im)


import  time
def draw(points):
    img = Image.fromarray(points)

    file_path2 = "{path}/ncData/{filename}".format(#注意网站和程序的目录不一样
       path=sys.path[0],
       filename= str(time.time())+'.tif'
    )
    img.save(file_path2)


def main(shapefile_path, raster_path):
    # Load the source data as a gdalnumeric array
    srcArray = gdalnumeric.LoadFile(raster_path)

    # Also load as a gdal image to get geotransform
    # (world file) info
    srcImage = gdal.Open(raster_path)
    geoTrans = srcImage.GetGeoTransform()

    # Create an OGR layer from a boundary shapefile
    shapef = ogr.Open(shapefile_path)
    lyr = shapef.GetLayer(os.path.split(os.path.splitext(shapefile_path)[0])[1])
    poly = lyr.GetNextFeature()

    # Convert the layer extent to image pixel coordinates
    minX, maxX, minY, maxY = lyr.GetExtent()
    ulX, ulY = world2Pixel(geoTrans, minX, maxY)
    lrX, lrY = world2Pixel(geoTrans, maxX, minY)

    # Calculate the pixel size of the new image
    pxWidth = int(lrX - ulX)
    pxHeight = int(lrY - ulY)

    clip = srcArray[ulY:lrY, ulX:lrX]

    #
    # EDIT: create pixel offset to pass to new image Projection info
    #
    xoffset = ulX
    yoffset = ulY
    print("Xoffset, Yoffset = ( %f, %f )" % (xoffset, yoffset))

    # Create a new geomatrix for the image
    geoTrans = list(geoTrans)
    geoTrans[0] = minX
    geoTrans[3] = maxY

    # Map points to pixels for drawing the
    # boundary on a blank 8-bit,
    # black and white, mask image.
    points = []
    pixels = []
    geom = poly.GetGeometryRef()
    pts = geom.GetGeometryRef(0)
    for p in range(pts.GetPointCount()):
        points.append((pts.GetX(p), pts.GetY(p)))
    for p in points:
        pixels.append(world2Pixel(geoTrans, p[0], p[1]))
    rasterPoly = Image.new("L", (pxWidth, pxHeight), 1)
    rasterize = ImageDraw.Draw(rasterPoly)
    rasterize.polygon(pixels, 0)
    mask = imageToArray(rasterPoly)

    # Clip the image using the mask
    # clip = gdalnumeric.choose(mask, \
    #                           (clip, 0)).astype(gdalnumeric.uint8)
    clip = gdalnumeric.choose(mask, (clip, 0))

    real = np.zeros(clip.shape)
    for i in range(len(clip)):
        for j in range(len(clip[i])):
            val = clip[i][j].real
            if (val < -4):
                val = -4
            real[i, j] = val

    dst_ds = gdal.GetDriverByName('GTiff').Create("hello86568.tif", 300, 300, 1, gdal.GDT_CFloat32)
    # dst_ds.SetGeoTransform([444720, 30, 0, 3751320, 0, -30])
    raster = np.zeros(real.shape, dtype=np.float32)
    dst_ds.GetRasterBand(1).WriteArray(real)
    # Once we're done, close properly the dataset
    dst_ds = None
    # This image has 3 bands so we stretch each one to make them
    # visually brighter
    # for i in range(3):
    clip2 = stretch(clip)

        # Save new tiff
    #
    #  EDIT: instead of SaveArray, let's break all the
    #  SaveArray steps out more explicity so
    #  we can overwrite the offset of the destination
    #  raster
    #
    ### the old way using SaveArray
    #
    # gdalnumeric.SaveArray(clip, "OUTPUT.tif", format="GTiff", prototype=raster_path)
    #
    ###
    #
    # gtiffDriver = gdal.GetDriverByName('GTiff')
    # if gtiffDriver is None:
    #     raise ValueError("Can't find GeoTiff Driver")
    # gtiffDriver.CreateCopy("beijing9.tif",
    #                        gdal_array.OpenArray(clip, prototype_ds=raster_path)
    #                        )
    # real = np.zeros(clip.shape)
    # for i in range(len(clip)):
    #     for j in range(len(clip[i])):
    #         val = clip[i][j].real
    #         real[i, j] = val
    # draw(real)

    # Save as an 8-bit jpeg for an easy, quick preview
    clip3 = clip2.astype(gdalnumeric.uint8)
    gdalnumeric.SaveArray(clip3, "beijing.jpg", format="JPEG")
    # misc.imsave("beijing7.png", clip2)

    gdal.ErrorReset()

import  sys
if __name__ == '__main__':
    raster = "D:\\FFOutput\\tsR1.tif"

    # 用于裁剪的shp
    # shp = "county"
    shp = "tsB.shp"
    # shapefile_path, raster_path
    # shapefile_path = 'temp/donghai.shp'
    # raster_path = 'temp/point1522361587.tif'
    main(shp, raster)
























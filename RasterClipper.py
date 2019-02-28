# -*- coding: UTF8 -*-

import operator
from osgeo import gdal, gdalnumeric, ogr, osr
from PIL import Image, ImageDraw
from functools import reduce

# 要裁剪的栅格图
# raster = "SatImage.tif"
raster = "D:\\FFOutput\\tsR1.tif"

# 用于裁剪的shp
# shp = "county"
shp = "D:\\FFOutput\\tsB"

# 裁剪后的栅格文件的名称
output = "clip"

# 此方法将栅格化的剪切shp转化为GDAL中使用的掩码
def imageToArray(i):
    """
    将PIL数组转换成gdalnumeric图像
    :param i:
    :return:
    """
    print(type(i))
    # a = gdalnumeric.fromstring(i.tostring(), 'b')
    a = gdalnumeric.fromstring(i.tobytes(), 'b')
    a.shape = i.im.size[1], i.im.size[0]
    return a


def arrayToImage(a):
    """
    将gdalnumeric数组转化为PIL图像
    :param a:
    :return:
    """
    # i = Image.fromstring('L', (a.shape[1], a.shape[0]), (a.astype('b')).tostring())
    print("-- arrayToImage'a -- > ", a)
    i = Image.frombytes('L', (a.shape[1], a.shape[0]), (a.astype('b')).tobytes())
    print("-- arrayToImage'i -- > ", i)
    return i


def world2Pixel(geoMatrix, x, y):
    """
    用gdal空间矩阵计算空间坐标点的像素位置
    :param geoMatrix:
    :param x:
    :param y:
    :return:
    """
    ulX = geoMatrix[0]
    ulY = geoMatrix[3]
    xDist = geoMatrix[1]
    yDist = geoMatrix[5]
    rtnX = geoMatrix[2]
    rtnY = geoMatrix[4]
    pixel = int((x - ulX) / xDist)
    line = int((ulY - y) / yDist)
    return (pixel, line)


def histogram(a, bins=range(0, 256)):
    """
    多维数组的直方图函数
    a = array
    bins = 要匹配的数字范围
    :param a:
    :param bins:
    :return:
    """
    fa = a.flat
    n = gdalnumeric.searchsorted(gdalnumeric.sort(fa), bins)
    n = gdalnumeric.concatenate([n, [len(fa)]])
    hist = n[1:] - n[:-1]
    return hist


def stretch(a):
    """
    在gdalnumeric数组图像上执行直方图拉伸
    :param a:
    :return:
    """
    hist = histogram(a)
    print("-- a --> ", a)
    im = arrayToImage(a)
    print("-- im --> ", im)
    lut = []
    for b in range(0, len(hist), 256):
        step = reduce(operator.add, hist[b:b+256]) / 255
        # 创建一个均衡查询表
        n = 0
        for i in range(256):
            lut.append(n / step)
            n = n + hist[i + b]
    im = im.point(lut)
    return imageToArray(im)


# 将源数据加载为一个gdalnumeric数组
srcArray = gdalnumeric.LoadFile(raster)

# 也可以加载为gdal图像以获取地理转换（世界文件）信息
srcImage = gdal.Open(raster)
geoTrans = srcImage.GetGeoTransform()
print("-- transformation -- > ", geoTrans)

# 从边界shp中创建一个OGR图层
shapef = ogr.Open("%s.shp" % shp)
# shapef = ogr.Open(shp)
lyr = shapef.GetLayer()
poly = lyr.GetNextFeature()

# 将图层范围转化为图像像素坐标
minX, maxX, minY, maxY = lyr.GetExtent()
print("-- minX --> ", minX)
print("-- maxX --> ", maxX)
print("-- minY --> ", minY)
print("-- maxY --> ", maxY)
ulX, ulY = world2Pixel(geoTrans, minX, maxY)
print("-- ulX --> ", ulX)
print("-- ulY --> ", ulY)
lrX, lrY = world2Pixel(geoTrans, maxX, minY)
print("-- lrX -->", lrX)
print("-- lrY -->", lrY)

# 计算新图像的像素大小
pxWidth = int(lrX - ulX)
pxHeight = int(ulY - lrY)
print("-- pxWidth --> ", pxWidth)
print("-- pxHeight --> ", pxHeight)

clip = srcArray[:, lrY:ulY, ulX:lrX]


# 为图像创建新的geomatrix
geoTrans = list(geoTrans)
geoTrans[0] = minX
geoTrans[3] = maxY


# 地图上的点到用在空白的8位黑白掩模图像上绘制边界的像素
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


# 用掩模裁剪图像
clip = gdalnumeric.choose(mask, (clip, 0)).astype(gdalnumeric.uint8)


# 这个图像有3个波段，因此我们得拉伸每一个波段使得他们有更好的可视化效果
for i in range(3):
    clip[i, :, :] = stretch(clip[i, :, :])

# 将ndvi另存为tiff
gdalnumeric.SaveArray(clip, "%s.tiff" % output, format="GTiff", prototype=raster)


# 将ndvi保存为8位的jpeg，以便于快速预览
clip = clip.astype(gdalnumeric.uint8)
gdalnumeric.SaveArray(clip, "%s.jpg" % output, format="JPEG")
























# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:38:07 2018

@author: MM
"""

import numpy as np
from osgeo import gdal

# Generate .asc header
def getHeader(ncols,nrows,xll,yll,cellsize,noData):
    header = "ncols        %s\n" % ncols
    header += "nrows        %s\n" % nrows
    header += "xllcorner    %s\n" % xll
    header += "yllcorner    %s\n" % yll
    header += "cellsize     %s\n" % cellsize
    header += "NODATA_value     %s" % noData
    return header

def openASC(filename,band=1):
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    return dsAsArray

def sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize,band=1):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur],xRes=cellsize,yRes=cellsize,resampleAlg=gdal.GRA_Mode)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%d",comments='')
    return dsAsArray

def averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize,band=1):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur],xRes=cellsize,yRes=cellsize,resampleAlg=gdal.GRA_Average)
    band = ds.GetRasterBand(band)
    dsAsArray = band.ReadAsArray()
    dsAsArray = np.clip(dsAsArray, 0, np.inf)
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%1.2f",comments='')
    return dsAsArray
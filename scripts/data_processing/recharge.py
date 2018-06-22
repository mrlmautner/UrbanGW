# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:33:24 2018

@author: MM
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from osgeo import gdal

xll = 457000
yll = 2108000
xur = 506000
yur = 2143000
cellsize = 500
ncols = int((506000-457000)/cellsize) # ncols
nrows = int((2143000-2108000)/cellsize) # nrows

def getHeader(ncols,nrows,xll,yll,cellsize,noData):
    header = "ncols        %s\n" % ncols
    header += "nrows        %s\n" % nrows
    header += "xllcorner    %s\n" % xll
    header += "yllcorner    %s\n" % yll
    header += "cellsize     %s\n" % cellsize
    header += "NODATA_value     %s" % noData
    return header

def openASC(filename):
    ds = gdal.Open(filename)
    band = ds.GetRasterBand(1)
    dsAsArray = band.ReadAsArray()
    return dsAsArray
    
LU = {}
RCH_Par = [0,0.5,1]
for year in [1984,1997,2003,2010,2012]:
    filename = r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\data\ModelInput\LU_' + str(year) + '.asc'
    LU_Array = openASC(filename)
    LU[str(year)] = LU_Array

header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
recharge = np.zeros((ncols,nrows))
for year in range(1984,2014):
    if year <= 1990:
        landusetxt = '1984'
    elif 1990 < year <= 2000:
        landusetxt = '1997'
    elif 2000 < year <= 2006:
        landusetxt = '2003'
    elif 2006 < year <= 2011:
        landusetxt = '2010'
    elif 2011 < year <= 2013:
        landusetxt = '2012'
        
    # Set effective precipitation percentages
    precipMult = np.zeros((LU[landusetxt].shape))
    precipMult[LU[landusetxt]==1] = RCH_Par[0]
    precipMult[LU[landusetxt]==2] = RCH_Par[1]
    precipMult[LU[landusetxt]==3] = RCH_Par[2]
    
    for month in range(1,13):
        
        filename = r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\data\ModelInput\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = openASC(filename)
        
        recharge = precipMult*precip/1000
        
        newfile = r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\data\ModelInput\RCH_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        np.savetxt(newfile, recharge, header=header, fmt="%1.5f",comments='')
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:33:24 2018

@author: MM
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal

os.chdir(r'E:\Tlalpan Model')

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500
ncols = int((xur-xll)/cellsize) # ncols
nrows = int((yur-yll)/cellsize) # nrows

def world2Pixel(gt, x, y):
  ulX = gt[0]
  ulY = gt[3]
  xDist = gt[1]
  yDist = gt[5]
  rtnX = gt[2]
  rtnY = gt[4]
  row = int((x - ulX) / xDist)
  col = int((ulY - y) / yDist)
  return (row, col)

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
 
#%%
LU = {}
RCH_Par = [1,0.5,1]
for year in [1984,1997,2003,2010,2012]:
    filename = r'Simulation_Wrapper\data\ModelInput\LU_' + str(year) + '.asc'
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
        
        filename = r'Simulation_Wrapper\data\ModelInput\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = openASC(filename)
        
        recharge = precipMult*precip/1000
        
        newfile = r'Simulation_Wrapper\data\ModelInput\75RCH_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        np.savetxt(newfile, recharge, header=header, fmt="%1.5f",comments='')
        
#%% Wells
WEL_Dict = {}

# Initialize dictionary with zero fluxes at layer 1, row 1, column 1
for period in range(0,360):
    WEL_Dict[period] = [[1,1,1,0]]

WEL_Array = np.loadtxt(r'Data\Pumping Wells\20180430_Pumping_AllDatasets_InModel_WatershedClip.csv', delimiter=',', skiprows=1, usecols=[1,2,3,4,5])
for w in range(0,WEL_Array.shape[0]):
    c = int(np.ceil((WEL_Array[w,0] - xll)/cellsize))
    r = int(np.ceil((yur - WEL_Array[w,1])/cellsize))
    for y in range(int(WEL_Array[w,3]-1984),int(WEL_Array[w,4]-1984)):
        for m in range(0,12):
            WEL_Dict[y*12+m].append([2,r,c,WEL_Array[w,2]])
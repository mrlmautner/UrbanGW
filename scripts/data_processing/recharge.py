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

os.chdir(r'F:\Tlalpan Model')

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
DAYS_M = [31,28,31,30,31,30,31,31,30,31,30,31]
LU = {}
RCH_Par = [1,0.5,1]
for year in [1984,1997,2003,2010,2012]:
    filename = r'Simulation_Wrapper\data\Input\LU_' + str(year) + '.asc'
    LU_Array = openASC(filename)
    LU[str(year)] = LU_Array
    
GEO = openASC(r'Simulation_Wrapper\data\Input\GEO_VM.asc')

#%% 1% Multiplier for lacustrine clays and 100% for all others
geoMult = 0.01*(GEO==1)+(GEO!=1)

#%%
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
for year in range(1984,2014):
    if year < 1991:
        landusetxt = '1984'
    elif 1991 <= year < 2001:
        landusetxt = '1997'
    elif 2001 <= year < 2005:
        landusetxt = '2003'
    elif 2005 <= year < 2011:
        landusetxt = '2010'
    elif 2011 <= year < 2014:
        landusetxt = '2012'
        
    # Set effective precipitation percentages
    precipMult1 = np.zeros((LU[landusetxt].shape))
    precipMult2 = np.zeros((LU[landusetxt].shape))
    precipMult3 = np.zeros((LU[landusetxt].shape))
    precipMult1[LU[landusetxt]==1] = 1
    precipMult2[LU[landusetxt]==2] = 1
    precipMult3[LU[landusetxt]==3] = 1
    
    for month in range(1,13):
        
        filename = r'Simulation_Wrapper\data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = openASC(filename)
        
        # Apply 1% recharge rate on the clay layer and separate by Recharge potential by LU type
        recharge1 = geoMult*precipMult1*precip/1000/DAYS_M[month-1]
        recharge2 = geoMult*precipMult2*precip/1000/DAYS_M[month-1]
        recharge3 = geoMult*precipMult3*precip/1000/DAYS_M[month-1]
        
        newfile1 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH1_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile2 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH2_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile3 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH3_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        np.savetxt(newfile1, recharge1, header=header, fmt="%1.5f",comments='')
        np.savetxt(newfile2, recharge2, header=header, fmt="%1.5f",comments='')
        np.savetxt(newfile3, recharge3, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for first time step
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    for month in range(1,13):
        filename = r'Simulation_Wrapper\data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = openASC(filename)
        avgPrecip += precip
        nprecip += 1

avgPrecip = avgPrecip/nprecip

precipMult1 = np.zeros((LU['1984'].shape))
precipMult2 = np.zeros((LU['1984'].shape))
precipMult3 = np.zeros((LU['1984'].shape))
precipMult1[LU['1984']==1] = 1
precipMult2[LU['1984']==2] = 1
precipMult3[LU['1984']==3] = 1

recharge1 = geoMult*precipMult1*avgPrecip/1000/31
recharge2 = geoMult*precipMult2*avgPrecip/1000/31
recharge3 = geoMult*precipMult3*avgPrecip/1000/31

newfile1 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH1_AVG.asc'
newfile2 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH2_AVG.asc'
newfile3 = r'Simulation_Wrapper\data\Input\RCH\ClayMult\RCH3_AVG.asc'

np.savetxt(newfile1, recharge1, header=header, fmt="%1.5f",comments='')
np.savetxt(newfile2, recharge2, header=header, fmt="%1.5f",comments='')
np.savetxt(newfile3, recharge3, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for monthly
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((12,nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    
    for month in range(1,13):
        filename = r'Simulation_Wrapper\data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = openASC(filename)
        avgPrecip[month-1] += precip

    nprecip += 1

for month in range(0,12):
    avgPrecip[month] = avgPrecip[month]/nprecip
    newfile1 = r'Simulation_Wrapper\data\Input\RCH\MTHLY\AVG_' + '{num:02d}'.format(num=month) + '.asc'
    np.savetxt(newfile1, avgPrecip[month], header=header, fmt="%1.5f",comments='')

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 16:38:07 2018

@author: MM
"""
import os
import numpy as np
from osgeo import gdal

os.chdir(r'C:\Users\MM\Google Drive\UrbanGW')
#os.chdir(r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper')

## Initialize model extent
# Lower left corner: (455000, 2107000)
# Upper right corner: (539000, 2175000)
xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500
ncols = int((xur-xll)/cellsize) # ncols
nrows = int((yur-yll)/cellsize) # nrows

# Generate .asc header
def getHeader(ncols,nrows,xll,yll,cellsize,noData):
    header = "ncols        %s\n" % ncols
    header += "nrows        %s\n" % nrows
    header += "xllcorner    %s\n" % xll
    header += "yllcorner    %s\n" % yll
    header += "cellsize     %s\n" % cellsize
    header += "NODATA_value     %s" % noData
    return header

def sample_Model_Griddata(filename,newfile,header):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur],xRes=cellsize,yRes=cellsize,resampleAlg=gdal.GRA_Mode)
    band = ds.GetRasterBand(1)
    dsAsArray = band.ReadAsArray()
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%d",comments='')
    return dsAsArray

def average_Model_Griddata(filename,newfile,header):
    ds = gdal.Warp(newfile,filename,outputBounds=[xll, yll, xur, yur],xRes=cellsize,yRes=cellsize,resampleAlg=gdal.GRA_Average)
    band = ds.GetRasterBand(1)
    dsAsArray = band.ReadAsArray()
    dsAsArray = np.clip(dsAsArray, 0, np.inf)
    
    np.savetxt(newfile, dsAsArray, header=header, fmt="%1.2f",comments='')
    return dsAsArray

#%%
# Precipitation rasters
for year in range(1984,2014):
    for month in range(1,13):
        filename = r'data\RawFiles\Precipitation\Precip_VM\Precip_mm_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile = r'data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
        dsAsArray = average_Model_Griddata(filename,newfile,header)

# Land Use rasters
for year in [1984,1997,2003,2010,2012,2016]:
    filename = r'data\RawFiles\LandUse\LU_' + str(year) + '_ZONE.asc'
    newfile = r'data\Input\LU_' + str(year) + '.asc'
    header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
    dsAsArray = sample_Model_Griddata(filename,newfile,header)
    
# Geology raster
filename = r'data\RawFiles\Geology\GEOLOGY_ZONES_M.asc'
newfile = r'data\Input\GEO_VM.asc'
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = sample_Model_Griddata(filename,newfile,header)

# DEM raster
filename = r'data\RawFiles\TopoRasters\DEM_Clipped.asc'
newfile = r'data\Input\DEM_VM.asc'
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = average_Model_Griddata(filename,newfile,header)

#%%

# Model 
filename = r'data\RawFiles\MODEL_ACTIVE.asc'
newfile = r'data\Input\ACTIVE.asc'
header = getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = sample_Model_Griddata(filename,newfile,header)
        
## Well data
#wellArray = np.loadtxt(r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\data\RawFiles\PumpingWells\PUMPING_M3_D_MC.csv', delimiter=',', skiprows=1, usecols=[3,4,7]) # (m3/day)
#PUMP_annual = np.zeros((30,3))
#PUMP_annual[:,0] = np.arange(1984,2014)
#for w in range(3003):
#    for i in range(0,int(wellArray[w,1] - wellArray[w,0])):
#        PUMP_annual[int(wellArray[w,0]+i-1984),1] += 1
#        PUMP_annual[int(wellArray[w,0]+i-1984),2] += wellArray[w,2]
#np.savetxt(r'C:\Users\MM\Google Drive\Tlalpan\Simulation_Wrapper\data\ModelInput\PUMP_annual.csv', PUMP_annual, delimiter=',')
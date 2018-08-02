# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 16:41:28 2018

@author: MM
"""
from gwscripts.dataprocessing import gengriddata as gen

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

#%%
# Precipitation rasters
for year in range(1984,2014):
    for month in range(1,13):
        filename = r'data\RawFiles\Precipitation\Precip_VM\Precip_mm_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile = r'data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
        dsAsArray = gen.averageGD(filename,newfile,header)

#%%
# Land Use rasters
for year in [1985,1990,1995,2000,2005,2010,2015]:
    filename = r'data_raw\LandUse\LU_' + str(year) + '_ZONE.asc'
    newfile = r'data\Input\LU_' + str(year) + '.asc'
    header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
    dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

#%%
# Geology raster
filename = r'data_raw\GEOLOGY_ZONES.asc'
newfile = r'data_output\GEO_VM.asc'
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

#%%
# DEM raster
filename = r'data_raw\DEM_Clipped.asc'
newfile = r'data_output\DEM_VM.asc'
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

#%%
# Active raster 
filename = r'data_raw\MODEL_ACTIVE.asc'
newfile = r'data_output\ACTIVE_VM.asc'
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)
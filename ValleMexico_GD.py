# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 16:41:28 2018

@author: MM
"""
from gwscripts.dataprocessing import gengriddata as gen
import numpy as np
import pandas as pd

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
#
##%%
## Precipitation rasters
#for year in range(1984,2014):
#    for month in range(1,13):
#        filename = r'data\RawFiles\Precipitation\Precip_VM\Precip_mm_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        newfile = r'data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#        dsAsArray = gen.averageGD(filename,newfile,header)
#
##%%
## Land Use rasters
#for year in [1985,1990,1995,2000,2005,2010,2015]:
#    filename = r'data_raw\LandUse\LU_' + str(year) + '_ZONE.asc'
#    newfile = r'data\Input\LU_' + str(year) + '.asc'
#    header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#    dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

#%%
# Geology raster
filename = r'data_raw\GEO_VM_INEGI-TR.asc'
newfile = r'data_output\GEO_VM.asc'
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

##%%
## Municipality raster
#filename = r'data_raw\MUN.asc'
#newfile = r'data_output\MUN_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)
#
##%%
## IH raster
#filename = r'data_raw\IH_OK_Combined_1984.asc'
#newfile = r'data_output\IH_1984.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

##%%
## DEM raster
#filename = r'data_raw\DEM_Clipped.asc'
#newfile = r'data_output\DEM_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

###%%
## THICKNESS raster
#filename = r'data_raw\THICK_LYR2.asc'
#newfile = r'data_output\THICK2_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)
#np.savetxt(newfile, dsAsArray, header=header, fmt="%1.2f",comments='')

##%%
## Active raster 
#filename = r'data_raw\MODEL_ACTIVE.asc'
#newfile = r'data_output\ACTIVE_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

##%%
## Well Data
#years = np.arange(1984,2014)
#dataset = years.copy()
#dataset[0:3] = 1987
#wellData = wellFeatures = pd.DataFrame(columns=['DATASET','X','Y','FLOW_M3-S','START_YR','END_YR','ID','S_PER','E_PER','S_DAY','E_DAY','Flow_m3','CVE_ENT','CVE_MUN','NOM_MUN','CODE'])
#for y in range(0,int(len(dataset))):
#    filename = r'data_raw\SACM_1987-2017\SACM_' + str(dataset[y]) + '.csv'
#    wellData = gen.generateWellData(filename, wellData, years[y])
#
#wellData.to_csv('data_output\wells\PUMP_S.csv', columns=['DATASET','X','Y','FLOW_M3-S','START_YR','END_YR','ID','S_PER','E_PER','S_DAY','E_DAY','Flow_m3','CVE_ENT','CVE_MUN','NOM_MUN','CODE'])
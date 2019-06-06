# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 16:41:28 2018

@author: MM
"""
from gwscripts.dataprocessing import gengriddata as gen
import os
from pathlib import Path

# Set current working directory to script folder
os.chdir(os.path.dirname(__file__))

# Set cwd to repository folder
path = Path(os.getcwd())
os.chdir(path)

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

## Precipitation rasters
#for year in range(1984,2014):
#    for month in range(1,13):
#        filename = r'data\RawFiles\Precipitation\Precip_VM\Precip_mm_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        newfile = r'data\Input\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#        dsAsArray = gen.averageGD(filename,newfile,header)

## Geology raster
#filename = Path.cwd()/'data_raw'/'GEO_VM_INEGI-TR.asc'
#newfile = Path.cwd()/'data_processed'/'GEO_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

## Municipality raster
#filename = r'data_raw\MUN.asc'
#newfile = r'data_processed\MUN_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

# IH raster
filename = str(Path.cwd()/'data_raw'/'IH-OK-Combo-DEM_1984_LT2750.asc')
newfile = str(Path.cwd()/'data_processed'/'IH_1984.asc')
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

## DEM raster
#filename = r'data_raw\DEM_Clipped.asc'
#newfile = r'data_processed\DEM_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

## THICKNESS raster
#filename = str(Path.cwd()/'data_raw'/'THICK2_VM_DRNS.asc')
#newfile = str(Path.cwd()/'data_processed'/'THICK2_VM.asc')
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.averageGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

## Active raster 
#filename = r'data_raw\MODEL_ACTIVE.asc'
#newfile = r'data_processed\ACTIVE_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#dsAsArray = gen.sampleGD(filename,newfile,header,xll,yll,xur,yur,cellsize)

## WETDRY raster
#LYR1FILE = r'data_processed\ACTIVE_VM_LYR1.asc'
#LYR2FILE = r'data_processed\ACTIVE_VM_LYR2.asc'
#newfile = r'data_processed\WETDRY_VM.asc'
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#lyr1 = np.loadtxt(LYR1FILE,skiprows=6)
#lyr2 = np.loadtxt(LYR2FILE,skiprows=6)
#wetdry = (lyr1 - lyr2) * 5
#np.savetxt(newfile, wetdry, header=header, fmt="%d",comments='')

## Well Data
#years = np.arange(1984,2014)
#dataset = years.copy()
#dataset[0:3] = 1987
#wellData = wellFeatures = pd.DataFrame(columns=['DATASET','X','Y','FLOW_M3-S','START_YR','END_YR','ID','S_PER','E_PER','S_DAY','E_DAY','Flow_m3','CVE_ENT','CVE_MUN','NOM_MUN','CODE'])
#for y in range(0,int(len(dataset))):
#    filename = r'data_raw\SACM_1987-2017\SACM_' + str(dataset[y]) + '.csv'
#    wellData = gen.generateWellData(filename, wellData, years[y])
#
#wellData.to_csv('data_processed\wells\PUMP_S.csv', columns=['DATASET','X','Y','FLOW_M3-S','START_YR','END_YR','ID','S_PER','E_PER','S_DAY','E_DAY','Flow_m3','CVE_ENT','CVE_MUN','NOM_MUN','CODE'])
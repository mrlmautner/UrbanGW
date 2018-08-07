# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 13:33:24 2018

@author: MM
"""
import numpy as np
from gwscripts.dataprocessing import gengriddata as gen

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500
ncols = int((xur-xll)/cellsize) # ncols
nrows = int((yur-yll)/cellsize) # nrows

rchFolder = 'UrbanGW//data_output//recharge//'
 
#%%
DAYS_M = [31,28,31,30,31,30,31,31,30,31,30,31]
LU = {}
for year in [1985,1990,1995,2000,2005,2010,2015]:
    filename = r'UrbanGW\data_output\LU_' + str(year) + '.asc'
    LU_Array = gen.openASC(filename)
    LU[str(year)] = LU_Array
    
GEO = gen.openASC(r'UrbanGW\data_output\GEO_VM.asc')

#%% 1% Multiplier for lacustrine clays and 100% for all others
geoMult = 0.01*(GEO==1)+(GEO!=1)

#%%
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
for year in range(1984,2014):
    if year < 1988:
        landusetxt = '1985'
    elif 1988 <= year < 1993:
        landusetxt = '1990'
    elif 1993 <= year < 1998:
        landusetxt = '1995'
    elif 1998 <= year < 2003:
        landusetxt = '2000'
    elif 2003 <= year < 2008:
        landusetxt = '2005'
    elif 2008 <= year < 2013:
        landusetxt = '2010'
    elif 2013 <= year < 2014:
        landusetxt = '2015'
    
    for month in range(1,13):
        
        filename = rchFolder + 'Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        
        # Apply 1% recharge rate on the clay layer and separate by Recharge potential by LU type
        recharge = geoMult*precip/1000/DAYS_M[month-1]
        
        newfile = rchFolder + 'ClayMult\\RCH_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        np.savetxt(newfile, recharge, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for first time step
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    for month in range(1,13):
        filename = rchFolder + 'Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        avgPrecip += precip
        nprecip += 1

avgPrecip = avgPrecip/nprecip

recharge = geoMult*avgPrecip/1000/30

newfile = rchFolder + 'ClayMult\\RCH_AVG.asc'

np.savetxt(newfile, recharge, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for monthly
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((12,nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    
    for month in range(1,13):
        filename = rchFolder + 'Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        avgPrecip[month-1] += precip

    nprecip += 1

for month in range(0,12):
    avgPrecip[month] = avgPrecip[month]/nprecip
    newfile1 = rchFolder + 'MTHLY\\AVG_' + '{num:02d}'.format(num=month) + '.asc'
    np.savetxt(newfile1, avgPrecip[month], header=header, fmt="%1.5f",comments='')

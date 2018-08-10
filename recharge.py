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
    
GEO = gen.openASC('UrbanGW//data_output//GEO_VM.asc',1)

#%% 1% Multiplier for lacustrine clays and 100% for all others
geoMult = 0.01*(GEO==1)+(GEO!=1)

#%%
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
for year in range(1984,2014):    
    for month in range(1,13):
        
        filename = rchFolder + 'Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        
        # Apply 1% recharge rate on the clay layer and separate by Recharge potential by LU type
        recharge = geoMult*precip/1000/DAYS_M[month-1]
        
        newfile = rchFolder + 'claymult\\RCH_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
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

newfile = rchFolder + 'claymult\\RCH_AVG.asc'

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
    newfile1 = rchFolder + 'mthly\\AVG_' + '{num:02d}'.format(num=month) + '.asc'
    np.savetxt(newfile1, avgPrecip[month], header=header, fmt="%1.5f",comments='')

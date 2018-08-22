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

rchFolder = 'data_output\\recharge\\'
 
#%%
DAYS_M = [31,28,31,30,31,30,31,31,30,31,30,31]
    
GEO = gen.openASC('data_output\GEO_VM.asc',1)

#%% 1% Multiplier for lacustrine clays and 100% for all others
geoMult = 0.01*(GEO==1)+(GEO!=1)

##%%
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#for year in range(1984,2014):    
#    for month in range(1,13):
#        
#        filename = rchFolder + 'precip\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        precip = gen.openASC(filename)
#        
#        # Apply 1% recharge rate on the clay layer and separate by Recharge potential by LU type
#        recharge = geoMult*precip/1000/DAYS_M[month-1]
#        
#        newfile = rchFolder + 'claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        np.savetxt(newfile, recharge, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for first time step
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
totPrecip = np.zeros((nrows,ncols))

for year in range(1984,2014):
    for month in range(1,13):
        filename = rchFolder + 'claymult\\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        totPrecip += precip*DAYS_M[month-1]

avgPrecip = totPrecip/30 # recharge divided by mm, 30 years, 365 days/year

newfile = rchFolder + 'claymult\\PrecipCM_AVGyr.asc'

np.savetxt(newfile, avgPrecip, header=header, fmt="%1.5f",comments='')

##%% Create Average Precipitation RCH file for monthly
#header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
#avgPrecip = np.zeros((12,nrows,ncols))
#nprecip = 0
#
#for year in range(1984,2014):
#    
#    for month in range(1,13):
#        filename = rchFolder + 'Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
#        precip = gen.openASC(filename)
#        avgPrecip[month-1] += precip
#
#    nprecip += 1
#
#for month in range(0,12):
#    avgPrecip[month] = avgPrecip[month]/nprecip
#    newfile1 = rchFolder + 'mthly\\AVG_' + '{num:02d}'.format(num=month) + '.asc'
#    np.savetxt(newfile1, avgPrecip[month], header=header, fmt="%1.5f",comments='')

#%% Determine precipitation and recharge for SS assumptions
Precip = gen.openASC(r'data_output\recharge\claymult\PrecipCM_AVG.asc')
Uper = gen.openASC(r'data_output\landuse\LU-1985-URBAN.asc')*Precip
Nper = gen.openASC(r'data_output\landuse\LU-1985-NATURAL.asc')*Precip
Wper = gen.openASC(r'data_output\landuse\LU-1985-WATER.asc')*Precip

Precip *= 500**2
Uper *= 0.01*500**2
Nper *= 0.30*500**2
Wper *= 0.5*500**2

P = np.sum(np.sum(Precip))
U = np.sum(np.sum(Uper))
N = np.sum(np.sum(Nper))
W = np.sum(np.sum(Wper))
R = U+N+W
ratio = R/P
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
 
#%%
DAYS_M = [31,28,31,30,31,30,31,31,30,31,30,31]
LU = {}
RCH_Par = [1,1,1]
for year in [1984,1997,2003,2010,2012]:
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
    
    # Set effective precipitation percentages
    precipMult1 = np.zeros((LU[landusetxt].shape))
    precipMult2 = np.zeros((LU[landusetxt].shape))
    precipMult3 = np.zeros((LU[landusetxt].shape))
    precipMult1[LU[landusetxt]==1] = 1
    precipMult2[LU[landusetxt]==2] = 1
    precipMult3[LU[landusetxt]==3] = 1
    
    for month in range(1,13):
        
        filename = r'UrbanGW\data_output\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        
        # Apply 1% recharge rate on the clay layer and separate by Recharge potential by LU type
        recharge1 = geoMult*precipMult1*precip/1000/DAYS_M[month-1]
        recharge2 = geoMult*precipMult2*precip/1000/DAYS_M[month-1]
        recharge3 = geoMult*precipMult3*precip/1000/DAYS_M[month-1]
        
        newfile1 = r'UrbanGW\data_output\RCH\ClayMult\RCH1_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile2 = r'UrbanGW\data_output\RCH\ClayMult\RCH2_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        newfile3 = r'UrbanGW\data_output\RCH\ClayMult\RCH3_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        np.savetxt(newfile1, recharge1, header=header, fmt="%1.5f",comments='')
        np.savetxt(newfile2, recharge2, header=header, fmt="%1.5f",comments='')
        np.savetxt(newfile3, recharge3, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for first time step
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    for month in range(1,13):
        filename = r'UrbanGW\data_output\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
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

newfile1 = r'UrbanGW\data_output\RCH\ClayMult\RCH1_AVG.asc'
newfile2 = r'UrbanGW\data_output\RCH\ClayMult\RCH2_AVG.asc'
newfile3 = r'UrbanGW\data_output\RCH\ClayMult\RCH3_AVG.asc'

np.savetxt(newfile1, recharge1, header=header, fmt="%1.5f",comments='')
np.savetxt(newfile2, recharge2, header=header, fmt="%1.5f",comments='')
np.savetxt(newfile3, recharge3, header=header, fmt="%1.5f",comments='')

#%% Create Average Precipitation RCH file for monthly
header = gen.getHeader(ncols,nrows,xll,yll,cellsize,-99999)
avgPrecip = np.zeros((12,nrows,ncols))
nprecip = 0

for year in range(1984,2014):
    
    for month in range(1,13):
        filename = r'UrbanGW\data_output\RCH\Precip_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        precip = gen.openASC(filename)
        avgPrecip[month-1] += precip

    nprecip += 1

for month in range(0,12):
    avgPrecip[month] = avgPrecip[month]/nprecip
    newfile1 = r'UrbanGW\data_output\RCH\MTHLY\AVG_' + '{num:02d}'.format(num=month) + '.asc'
    np.savetxt(newfile1, avgPrecip[month], header=header, fmt="%1.5f",comments='')

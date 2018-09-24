# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import os
import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import calendar
import time
from gwscripts.dataprocessing import gengriddata as gen
from gwscripts.flopymodel import flomodelvm as mod

timestart = time.time()
print('Processing data...')

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

STRT_YEAR = 1984
END_YEAR = 2014

# Load datasets
ACTIVE = gen.openASC('data_output\ACTIVE_VM.asc')
GEO = gen.openASC('data_output\GEO_VM.asc')
DEM = gen.openASC('data_output\DEM_VM.asc')
IH = gen.openASC('data_output\IH_1984.asc')

# Assign name and create modflow model object
modelname = 'VM_Test'

#%% Process parameter data
# COL1: Phase period starting stress period
PHASE_PER = [2, 121, 241, 361]
# COL2 : Phase year
LU_PAR = [1985, 1990, 2000, 2010]
# COL3 : Phase well pumping multiplier 
WEL_PAR = [0.58, 1.33,1.446,1.687]
#COL3 : Phase distribution system leak multiplier
LEAK_PAR = np.array([40.94, 40.94, 40.94, 40.94])
# COL4 : Phase LID increase multiplier
LID_PAR = [1, 1, 1, 1]
# COL5 :  
PhaseParams = np.array([PHASE_PER,LU_PAR,WEL_PAR,LEAK_PAR,LID_PAR])

# COL1 : Zone hydraulic conductivity vertical anisotropy
VK_PAR = [100,100,10,0.1,0.01,10]
# COL2 : Zone specific storage
SS_PAR = [6.56E-02,3.28E-04,3.28E-04,1.64E-05,1.64E-05,3.28E-06]
# COL3 : Zone hydraulic conductivity
HK_PAR = [4.32E-03,43.20,0.432,4.32,0.0432,8.64E-05] # m/d
# COL4 : Zone specific yield
SY_PAR = [0.06,0.15,0.15,0.30,0.01,0.01]
ZoneParams = np.array([HK_PAR,VK_PAR,SS_PAR,SY_PAR])

RCH_PAR = [1.00E-02, 0.7051, 50.00E-02] # Recharge multiplier for urban, natural, and water cover

ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                 STRT_YEAR,END_YEAR,ACTIVE,GEO,DEM,IH,
                                                 ZoneParams)

print('Basic, Discretization, and Layer packages generated in',str(time.time()-timestart),'seconds')

#%% Assign Land Use Types
#LU = {'0':{},'1':{},'2':{}}
#LU_PAR = [1985, 1990, 2000, 2010]
#for n in range(0,4):
#    for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
#        filename = r'C:\Users\MM\Google Drive\UrbanGW\data_output\landuse\LU-' + str(LU_PAR[n]) + '-' + LUtype + '.asc'
#        LU[str(l)][str(LU_PAR[n])] = gen.openASC(filename)


#%% Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
newtime = time.time()

PUMP_PARAM = [1,WEL_PAR[0],WEL_PAR[1],WEL_PAR[2],WEL_PAR[0],WEL_PAR[1],WEL_PAR[2],WEL_PAR[3]]
start = [1984-1/12,1984-1/12,1984,1995,1984-1/12,1984,1995,2005]
end = [1984,1984,1995,2005,1984,1995,2005,2014]

PUMP_array = np.loadtxt(r'data_output\wells\PUMP_CS.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
WEL_DICT = mod.addNewWells(PUMP_array,LYR=1)

# Add supply wells
for i, wellset in enumerate(['PUMP_C1984','PUMP_S05_Q','PUMP_S05_Q','PUMP_S05_Q','PUMP_RC_Q','PUMP_RC_Q','PUMP_RC_Q','PUMP_RC_Q']):
    PUMP_array = np.loadtxt(r'data_output\wells\\' + wellset + '.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
    
    WEL_DICT = mod.addNewWells(New_WEL=PUMP_array,LYR=1,WEL_Dict=WEL_DICT,WEL_mult=PUMP_PARAM[i],S_YR=start[i],E_YR=end[i])

# Add leak wells
start = [1984-1/12,1984,1995,2005]
end = [1984,1995,2005,2014]
for i, leakset in enumerate(['1985','1990','2000','2010']):
    LEAK_array = np.loadtxt(r'data_output\leak\LEAK_' + leakset + '.csv',
                                              delimiter=',', skiprows=1, usecols=[2,1,7]) # pumping in m3 per day
    LEAK_array = np.insert(LEAK_array, 2, start[i], axis=1)
    LEAK_array = np.insert(LEAK_array, 3, end[i], axis=1)
    WEL_DICT = mod.addNewWells(LEAK_array,LYR=1,WEL_Dict=WEL_DICT,WEL_mult=LEAK_PAR[i],coordType='rc')

## Add the difference between the installed WWTP minus the actual treatment quantity (m3/s)
#WWTP_array = np.loadtxt(r'data_output\scenarios\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6])
#WWTP_array[:,3] = WWTP_array[:,2] - WWTP_array[:,3]
#WWTP_array = np.insert(WWTP_array, 2, 1984, axis=1)
#WWTP_array[WWTP_array[:,3]<0.1,2] = 1984+5
#WWTP_array[WWTP_array[:,3]>=0.1,2] = 1984+10
#WWTP_array[:,3] = np.ones(WWTP_array.shape[0])*2014
#WEL_DICT = mod.addNewWells(WWTP_array,LYR=1,WEL_Dict=WEL_DICT,WEL_mult=60*60*24) # Mult used to convert to m3/d

#%% Recharge Basins
#Basin_Array = np.loadtxt(r'data\Input\Scenarios\RechargeBasins_Elevation.csv', delimiter=',', skiprows=1, usecols=[13,15,16])
#bLoc = [[22,121],[48,122],[54,92],[66,125],[84,120]]
#for b in range(0,5):
##    randBasin = np.random.randint(0,Basin_Array.shape[0])
##    c = int(np.floor((Basin_Array[randBasin,1] - xll)/cellsize))
##    r = int(np.floor((yur - Basin_Array[randBasin,2])/cellsize))
#    r = bLoc[b][0]
#    c = bLoc[b][1]
#    
#    for per in range(5*12,360):
#        # Add injection equal to treatment capacity for each parameter period
#        WEL_Dict[per].append([0,r,c,0.9648*86400]) # 1 m3/s recharge basins (35 cfs)

wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)

print('WEL_Dict generated in',str(time.time()-newtime),'seconds')

#%% Recharge
newtime = time.time()

start = [0,1984,1995,2005]
end = [0,1995,2005,2014]
rch = {0:[0.01,0.267,0.5],1:RCH_PAR,2:RCH_PAR,3:RCH_PAR}

RCH_DICT = {}
Precip_Dict = {}

filename = r'data_output\recharge\claymult\PrecipCM_AVG.asc'
Precip_Dict[0] = gen.openASC(filename)

for year in range(int(STRT_YEAR),int(END_YEAR)):
    for month in range(1,13):
        per = (year-1984)*12+month
    
        filename = r'data_output\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        Precip_Dict[per] = gen.openASC(filename)
        
for i, LUset in enumerate(['1985','1990','2000','2010']):
    LU = [0]*3
    
    for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
        filename = r'data_output\landuse\LU-' + LUset + '-' + LUtype + '.asc'
        LU[l] = gen.openASC(filename)
    
    RCH_DICT = mod.addRecharge(LU_arrays=LU,PRECIP=Precip_Dict,S_YR=start[i],E_YR=end[i],RCH_Dict=RCH_DICT,RCH_mult=rch[i])

rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)

print('RCH_Dict generated in',str(time.time()-newtime),'seconds')

oc, pcg = mod.outputControl(mf)
#oc = flopy.modflow.ModflowOc.load('ValleMexicoTRC.oc', mf)

# Run Model and post processing
# Write the MODFLOW model input files

print('Data processed in',str(time.time()-timestart),'seconds')

newtime = time.time()
print('Writing input file...')

mf.write_input()
print('Input file written in',str(time.time()-newtime),'seconds')

# Run the MODFLOW model
success, buff = mf.run_model()
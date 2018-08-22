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

H_INIT = bf.HeadFile('data_output\IH-from-SS.hds')

# Assign name and create modflow model object
modelname = 'VM_Test'

#%% Process parameter data
# COL1: Phase period starting stress period
PHASE_PER = [2, 121, 241, 361]
# COL2 : Phase year
LU_PAR = [1985, 1990, 2000, 2010]
# COL3 : Phase well pumping multiplier 
WEL_PAR = [0.64, 1.42, 1.49, 2.00]
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

RCH_PAR = [1.00E-02, 25.00E-02, 50.00E-02] # Recharge multiplier for urban, natural, and water cover

ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                 STRT_YEAR,END_YEAR,ACTIVE,GEO,DEM,H_INIT,
                                                 PhaseParams,ZoneParams)
#%% Assign Land Use Types
LU = {'0':{},'1':{},'2':{}}
LU_PAR = [1985, 1990, 2000, 2010]
for n in range(0,4):
    for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
        filename = r'C:\Users\MM\Google Drive\UrbanGW\data_output\landuse\LU-' + str(LU_PAR[n]) + '-' + LUtype + '.asc'
        LU[str(l)][str(LU_PAR[n])] = gen.openASC(filename)

##%% Wells
#WEL_Dict = {}
#WEL_map = np.zeros((136,168))
#
## Initialize dictionary with zero fluxes at layer 1, row 1, column 1
#for period in range(0,360):
#    WEL_Dict[period] = [[0,1,1,0]]
#
## Add pumping wells for all years in the study period (extend datasets backwards)
#WEL_Array = np.loadtxt(r'data\RawFiles\PumpingWells\20180430_Pumping_AllDatasets_InModel_WatershedClip_Abridged.csv',
#                                              delimiter=',', skiprows=1, usecols=[1,2,6,9,10]) # pumping in m3 per day
#for w in range(0,WEL_Array.shape[0]):
#    c = int(np.floor((WEL_Array[w,0] - xll)/cellsize))
#    r = int(np.floor((yur - WEL_Array[w,1])/cellsize))
#    
#    for per in range(int(WEL_Array[w,3]),int(WEL_Array[w,4])):
#        # Assign pumping multiplier for each parameter period
#        if per < PAR_LIM[0]:
#            param = 0
#        else:
#            if PAR_LIM[0] <= per < PAR_LIM[1]:
#                param = 0
#            elif PAR_LIM[1] <= per < PAR_LIM[2]:
#                param = 1
#            elif PAR_LIM[2] <= per < PAR_LIM[3]:
#                param = 2
#            elif PAR_LIM[3] <= per < PAR_LIM[4]:
#                param = 3
#            elif PAR_LIM[4] <= per:
#                param = 4
#            
#            WEL_Dict[per].append([1,r,c,WEL_Array[w,2]*WEL_PAR[param]])
#
## Add leak wells for all years in the study period
#for param in range(0,4):
#    for n in range(0,URBAN[param].shape[0]):
#        for per in range(PAR_LIM[param],PAR_LIM[param+1]):
#            WEL_Dict[per].append([1,URBAN[param][n,2],URBAN[param][n,3],LEAK_PAR[param]/URBAN[param].shape[0]])
#
#wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_Dict)
#
##%% Recharge
#RCH_Dict = {}
#precip = []
#
#for year in range(1984,2014):
#    for month in range(0,12):
#        
#        per = (year-1984)*12+month
#        RCH_Array = np.zeros((nrow,ncol))
#        
#        for lu in range(0,3):
#            filename = r'data\Input\RCH\ClayMult\RCH' + str(lu+1) + '_' + str(year) + '_' + '{num:02d}'.format(num=month+1) + '.asc'
#            precip = gen.openASC(filename)
#            RCH_Array += precip*RCH_PAR[lu]
#        
#        RCH_Dict[per] = RCH_Array
#        
#for lu in range(0,3):
#    filename = r'data\Input\RCH\RCH' + str(lu+1) + '_AVG.asc'
#    precip = gen.openASC(filename)
#    RCH_Array += precip*RCH_PAR[lu]
#    
#    RCH_Dict[0] = RCH_Array
#    
#rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_Dict)
#
##%% Output Control and Solver
## Add OC package to the MODFLOW model
#spd = {}
#data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
#spd[0,30] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference' ]
#for m in range(0,12):
#    for d in range(0,calendar.monthrange(1984,m+1)[1]):
#        spd[m,d] = data2record.copy()
#for y in range(0,30):
#    for m in range(1,12):
#        spd[y*12+m,(calendar.monthrange(1984+y,m+1)[1]-1)] = data2record.copy()
#for p in [6,20,32,41,54,78,90,102,114,128,138,150,162,175,187,198,213,225,235]:
#    spd[p,0] = data2record.copy()
#oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
#
## Add PCG package to the MODFLOW model
#pcg = flopy.modflow.ModflowPcg(mf)
#
##%% Run Model and post processing
## Write the MODFLOW model input files
#print('Data processed in',str(time.time()-timestart),'seconds')
#print('Writing input file...')
#mf.write_input()
#print('Input file written in',str(time.time()-timestart),'seconds')
#
## Run the MODFLOW model
#success, buff = mf.run_model()
#
#hds = bf.HeadFile(modelname+'.hds')
#h = hds.get_data(totim=1.0)
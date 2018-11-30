# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import time
from gwscripts.dataprocessing import gengriddata as gen
from gwscripts.flopymodel import flomodelvm as mod
import pickle

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
ACTIVE_LYR1 = gen.openASC('data_output\ACTIVE_VM_LYR1.asc')
ACTIVE_LYR2 = gen.openASC('data_output\ACTIVE_VM_LYR2.asc')
TH1 = gen.openASC('data_output\THICK1_VM.asc')
TH2 = gen.openASC('data_output\THICK2_VM.asc')
GEO = gen.openASC('data_output\GEO_VM.asc')
DEM = gen.openASC('data_output\DEM_VM.asc')
IH = gen.openASC('data_output\IH_1984.asc')
MUN = gen.openASC('data_output\MUN_VM.asc')

scenario = 'WWTP' # Historical' # Historical' # Test' # Leak' # Basin' # 
fixleak = 1
# Assign name and create modflow model object
modelname = 'model_output\VM_'+scenario

#%% Process parameter data
# COL1: Phase period starting stress period
PHASE_PER = [0, 132, 252, 360]
# COL2 : Phase year
LU_PAR = [1990, 2000, 2010]
# COL3 : Phase well pumping multiplier 
WEL_PAR = np.array([4.000E+00,3.903E+00,3.325E+00])
#COL3 : Phase distribution system leak multiplier
LEAK_PAR = np.array([1,1,1])
# COL4 : Phase LID increase multiplier
LID_PAR = [1, 1, 1]
# COL5 :  
PhaseParams = np.array([PHASE_PER,LU_PAR,WEL_PAR,LEAK_PAR,LID_PAR])

# COL1 : Zone hydraulic conductivity vertical anisotropy
VK_PAR = [100,100,10,1,0.1]
# COL2 : Zone specific storage
SS_PAR = [6.562E-02,1.960E-03,5.888E-03,6.496E-07,3.045E-07]
# COL3 : Zone hydraulic conductivity
HK_PAR = [4.320E-04,4.375E+02,6.995E-01,3.458E-01,8.640E-02] # m/d
# COL4 : Zone specific yield
SY_PAR = [0.06,0.15,0.15,0.30,0.01]
ZoneParams = np.array([HK_PAR,VK_PAR,SS_PAR,SY_PAR])

municipalities = np.unique(MUN)[1:]
#municipalities = dict.fromkeys(municipalities,{'UAREA':{'1990':{},'2000':{},'2010':{}},'OBJECTIVES':np.zeros(3)})

RCH_PAR = [1.00E-02, 1.215E-02, 50.00E-02] # Recharge multiplier for urban, natural, and water cover

ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                 STRT_YEAR,END_YEAR,ACTIVE_LYR1,ACTIVE_LYR2,GEO,DEM,IH,
                                                 ZoneParams,THICK1=TH1,THICK2=TH2)

print('Basic, Discretization, and Layer packages generated in',str(time.time()-timestart),'seconds')

#%% Land Use Type
# Fill a land use dictionary with the ARRAYs that represent the % of each land use cover in each cell
# and the LISTs that contain all the cells and percentages of each land use type
LU = {'1990':{},'2000':{},'2010':{}}     
for i, LUset in enumerate(['1990','2000','2010']):
    LU[LUset] = {'ARRAY':{},'LIST':{}}
    
    for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
        filename = r'data_output\landuse\LU-' + LUset + '-' + LUtype + '.asc'
        perarea =  gen.openASC(filename)
        LU[LUset]['ARRAY'][LUtype] = perarea
        
        LU[LUset]['LIST'][LUtype] = np.zeros((perarea.shape[0]*perarea.shape[1],5))
        
        l = 0
        for row in range(0,perarea.shape[0]):
            for col in range(0,perarea.shape[1]):
                if perarea[row,col] > 0.001:
                    LU[LUset]['LIST'][LUtype][l,2] = perarea[row,col]
                LU[LUset]['LIST'][LUtype][l,0] = col
                LU[LUset]['LIST'][LUtype][l,1] = row
                LU[LUset]['LIST'][LUtype][l,3] = 1-ACTIVE_LYR1[row,col] # 0 if clay layer, 1 if no clay layer
                LU[LUset]['LIST'][LUtype][l,4] = MUN[row,col]
                l += 1
        LU[LUset]['LIST'][LUtype] = LU[LUset]['LIST'][LUtype][LU[LUset]['LIST'][LUtype][:,2]>0,:]

winfofile = 'model_output\objective_data\LU_'+scenario+'.pickle'
with open(winfofile, 'wb') as handle:
    pickle.dump(LU, handle, protocol=pickle.HIGHEST_PROTOCOL)

#%% Recharge
newtime = time.time()

start = [1984,1995,2005]
end = [1995,2005,2014]

RCH_DICT = {}
Precip_Dict = {}

#filename = r'data_output\recharge\claymult\PrecipCM_AVG.asc'
#Precip_Dict[0] = gen.openASC(filename)

for year in range(int(STRT_YEAR),int(END_YEAR)):
    for month in range(1,13):
        per = (year-STRT_YEAR)*12+month-1
    
        filename = r'data_output\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        Precip_Dict[per] = gen.openASC(filename)

for i, LUset in enumerate(['1990','2000','2010']):
    RCH_DICT = mod.addRecharge(LU_arrays=LU[LUset]['ARRAY'],PRECIP=Precip_Dict,S_YR=start[i],E_YR=end[i],RCH_Dict=RCH_DICT,RCH_mult=RCH_PAR)

rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)

print('RCH_Dict generated in',str(time.time()-newtime),'seconds')

#%% Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
newtime = time.time()

WEL_DICT = {}
WEL_INFO = {}

# Add leak wells
mun = municipalities.copy()
mun = np.delete(mun,29)
mun = np.delete(mun,21)

start = [1984,1995,2005]
end = [1995,2005,2014]
LEAK_arrays = {'1990':{},'2000':{},'2010':{}}
for i, leakset in enumerate(['1990','2000','2010']):
    LEAK_MUN = np.loadtxt('data_output\leak\LEAK_TOT_MUN.csv',delimiter=',',skiprows=1,usecols=[0,i+1,4]) # Total leak per municipality m3/d
    LEAK_arrays[leakset] = np.zeros((LU[leakset]['LIST']['URBAN'].shape[0],3))
    j = 0
    
    for n,m in enumerate(mun):
        # Calculate the total urban area in normalized cell area (1 cell = 1) per municipality
        tempLeak = LU[leakset]['LIST']['URBAN'][(LU[leakset]['LIST']['URBAN'][:,4]==m),:4]
        u_area = tempLeak[:,2].sum()
        
        LperArea = LEAK_MUN[n,1]/u_area
        
        tempLeak[:,2] *= LperArea
        tempLeak[tempLeak[:,3]==1,2] *= 0.1 # apply 90% returns to sewer under clay layer
        LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),:] = tempLeak[:,:3]
        j += tempLeak.shape[0]
    
    LEAK_arrays[leakset] = LEAK_arrays[leakset][(LEAK_arrays[leakset][:,2]>5),:]
        
    LEAK_arrays[leakset] = np.insert(LEAK_arrays[leakset], 2, start[i], axis=1)
    LEAK_arrays[leakset] = np.insert(LEAK_arrays[leakset], 3, end[i], axis=1)
    
    WEL_DICT, WEL_INFO = mod.addNewWells(LEAK_arrays[leakset],LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=fixleak*LEAK_PAR[i],coordType='rc0',mun=MUN,wellType=1)

# Add supply wells, includes the fixleak and LEAK_MUN terms to reduce pumping when leaks are fixed
PUMP_PARAM = [WEL_PAR[0],WEL_PAR[1],WEL_PAR[2]]
start = [1984,1995,2005]
end = [1995,2005,2014]

PUMP_array = np.loadtxt(r'data_output\wells\PUMP_C.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,dateType='per',mun=MUN,munleak=LEAK_MUN,F=fixleak)
PUMP_array = np.loadtxt(r'data_output\wells\PUMP_S.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,dateType='per',mun=MUN,munleak=LEAK_MUN,F=fixleak)
 
for i, wellset in enumerate(['PUMP_RC_Q','PUMP_RC_Q','PUMP_RC_Q']):
    PUMP_array = np.loadtxt(r'data_output\wells\\' + wellset + '.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
    
    WEL_DICT, WEL_INFO = mod.addNewWells(New_WEL=PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=PUMP_PARAM[i],S_YR=start[i],E_YR=end[i],mun=MUN,munleak=LEAK_MUN,F=fixleak)

# Add the difference between the installed WWTP minus the actual treatment quantity (m3/s)
WWTP_array = np.loadtxt(r'data_output\scenarios\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6])
WWTP_array[:,3] = WWTP_array[:,2] - WWTP_array[:,3]
WWTP_array = np.insert(WWTP_array, 2, 1984, axis=1)
#WWTP_array[WWTP_array[:,3]<0.1,2] = 1984+5
#WWTP_array[WWTP_array[:,3]>=0.1,2] = 1984+10
WWTP_array[:,3] = np.ones(WWTP_array.shape[0])*2014
WEL_DICT, WEL_INFO = mod.addNewWells(WWTP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=60*60*24,mun=MUN,wellType=-1) # Mult used to convert to m3/d

##%% Recharge Basins
##Basin_Array = np.loadtxt(r'data\Input\Scenarios\RechargeBasins_Elevation.csv', delimiter=',', skiprows=1, usecols=[13,15,16])
#bLoc = [[22,121],[48,122],[54,92],[66,125],[84,120]]
#for b in range(0,5):
##    randBasin = np.random.randint(0,Basin_Array.shape[0])
##    c = int(np.floor((Basin_Array[randBasin,1] - xll)/cellsize))
##    r = int(np.floor((yur - Basin_Array[randBasin,2])/cellsize))
#    r = bLoc[b][0]
#    c = bLoc[b][1]
#    
#    for per in range(0,360):
#        # Add injection equal to treatment capacity for each parameter period
#        WEL_DICT[per].append([1,r,c,0.9648*86400]) # 1 m3/s recharge basins (35 cfs)
#        WEL_INFO[per].append([1,r,c,0.9648*86400,MUN[r,c],1])

wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)

print('WEL_Dict generated in',str(time.time()-newtime),'seconds')
newtime = time.time()

winfofile = 'model_output\objective_data\WEL_INFO_'+scenario+'.pickle'
with open(winfofile, 'wb') as handle:
    pickle.dump(WEL_INFO, handle, protocol=pickle.HIGHEST_PROTOCOL)

print('WEL_Dict saved in',str(time.time()-newtime),'seconds')    

#%%
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
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 16:36:38 2018

@author: MM
"""
import flopy
import numpy as np
import time
from gwscripts.dataprocessing import gengriddata as gen
from gwscripts.flopymodel import flomodelvm as mod

#def runOTCmodel(xll,yll):
timestart = time.time()
print('Processing data...')

xll = 250
yll = 250
xur = 5250
yur = 5250
cellsize = 500
ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns

STRT_YEAR = 1984
END_YEAR = 1987

# Load datasets
ACTIVE1 = gen.openASC('data_output\objTest\ACTIVE_VM.asc')
ACTIVE2 = np.ones((nrow,ncol))
GEO = gen.openASC('data_output\objTest\GEO_VM.asc')
DEM = gen.openASC('data_output\objTest\DEM_VM.asc')
IH_SS = DEM #gen.openASC('data_output\objTest\IH_1984.asc')
MUN = gen.openASC('data_output\objTest\MUN_VM.asc')

# Assign name and create modflow model object
modelname = 'OBJ_Test'

#%% Process parameter data
# COL2 : Phase year
LU_PAR = 1990
# COL3 : Phase well pumping multiplier 
WEL_PAR = 440
#COL3 : Phase distribution system leak multiplier
LEAK_PAR = 1330
# COL4 : Phase LID increase multiplier
LID_PAR = 1

# COL1 : Zone hydraulic conductivity vertical anisotropy
VK_PAR = [100,100,10,1,1]
# COL2 : Zone specific storage
SS_PAR = [6.56E-02,3.28E-04,3.28E-04,1.64E-05,1.64E-05]
# COL3 : Zone hydraulic conductivity
HK_PAR = [4.32E-03,43.20,0.432,4.32,0.0432] # m/d
# COL4 : Zone specific yield
SY_PAR = [0.06,0.15,0.15,0.30,0.01]
ZoneParams = np.array([HK_PAR,VK_PAR,SS_PAR,SY_PAR])

RCH_PAR = [0.003*42, 0.30*42, 0.50*42] # Recharge multiplier for urban, natural, and water cover

ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                             STRT_YEAR,END_YEAR,ACTIVE1,ACTIVE2,
                                                             GEO,DEM,IH_SS,ZoneParams,THICK1=50,THICK2=DEM-1700)

print('Basic, Discretization, and Layer packages generated in',str(time.time()-timestart),'seconds')

#%% Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
newtime = time.time()

PUMP_array = np.loadtxt(r'data_output\objTest\wells\PUMP_objTC.csv',
                                              delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
WEL_DICT,WEL_INFO = mod.addNewWells(New_WEL=PUMP_array,LYR=1,WEL_mult=WEL_PAR,dateType='per',coordType='rc',mun=MUN)

# Add leak wells
LEAK_array = np.loadtxt(r'data_output\objTest\leak\LEAK_' + str(LU_PAR) + '.csv',
                                          delimiter=',', skiprows=1, usecols=[2,1,7]) # pumping in m3 per day
LEAK_array = np.insert(LEAK_array, 2, STRT_YEAR, axis=1)
LEAK_array = np.insert(LEAK_array, 3, END_YEAR, axis=1)
WEL_DICT, WEL_INFO = mod.addNewWells(LEAK_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=LEAK_PAR,coordType='rc',mun=MUN,wellType=1)

# Add the difference between the installed WWTP minus the actual treatment quantity (m3/s)
WWTP_array = np.loadtxt(r'data_output\objTest\scenarios\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6])
WWTP_array[:,3] = WWTP_array[:,2] - WWTP_array[:,3]
WWTP_array = np.insert(WWTP_array, 2, 1984, axis=1)
WWTP_array[:,3] = np.ones(WWTP_array.shape[0])*1987
WEL_DICT, WEL_INFO = mod.addNewWells(WWTP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=60*60*24,mun=MUN,xll=xll,yur=yur,wellType=-1) # Mult used to convert to m3/d

#%% Recharge Basins
bLoc = [[5,5],[7,7]]
for b in range(0,2):
#    randBasin = np.random.randint(0,Basin_Array.shape[0])
#    c = int(np.floor((Basin_Array[randBasin,1] - xll)/cellsize))
#    r = int(np.floor((yur - Basin_Array[randBasin,2])/cellsize))
    r = bLoc[b][0]
    c = bLoc[b][1]
    
    for per in range(0,36):
        # Add injection equal to treatment capacity for each parameter period
        WEL_DICT[per].append([0,r,c,0.9648*86400]) # 1 m3/s recharge basins (35 cfs)
        WEL_INFO[per].append([0,r,c,0.9648*86400,MUN[r,c],1])
wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT, options='AUX IFACE')

print('WEL_Dict generated in',str(time.time()-newtime),'seconds')

#%% Recharge
newtime = time.time()

RCH_DICT = {}
Precip_Dict = {}

filename = r'data_output\objTest\recharge\claymult\PrecipCM_AVG.asc'
Precip_Dict[0] = gen.openASC(filename)

for year in range(int(STRT_YEAR),int(END_YEAR)):
    for month in range(1,13):
        per = (year-STRT_YEAR)*12+month
    
        filename = r'data_output\objTest\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
        Precip_Dict[per] = gen.openASC(filename)
        
LU = {}

for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
    filename = r'data_output\objTest\landuse\LU-' + str(LU_PAR) + '-' + LUtype + '.asc'
    LU[LUtype] = gen.openASC(filename)

RCH_DICT = mod.addRecharge(LU_arrays=LU,PRECIP=Precip_Dict,S_YR=STRT_YEAR,E_YR=END_YEAR,RCH_Dict=RCH_DICT,RCH_mult=RCH_PAR)

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
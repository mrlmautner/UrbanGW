# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
import calendar
import time

timestart = time.time()
print('Processing data...')

def initializeFM(modelname,xll,yll,xur,yur,cellsize,STRT_YEAR,END_YEAR,ACTIVE,GEO,DEM,ZoneParams):
    # modelname to set the file root 

    #mf = flopy.modflow.Modflow.load(r'ValleyMexico_zones.nam',exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')
    mf = flopy.modflow.Modflow(modelname, exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')
    
    # Model domain and grid definition
    ztop = DEM # Model top elevation (Digital Elevation Model)
    L1botm = ztop - 50 # Layer 1 bottom elevation
    L2botm = L1botm - 350 # Layer 2 bottom elevation
    L3botm = L2botm - 1500 # Layer 3 bottom elevation
    botm = [L1botm,L2botm,L3botm] # Model bottom elevation
    nlay = 3 # Number of layers
    ncol = int((xur-xll)/cellsize) # Number of rows
    nrow = int((yur-yll)/cellsize) # Number of columns
    delr = cellsize # Row height
    delc = cellsize # Column height
    
    # Time discretization
    nper = (END_YEAR - STRT_YEAR)*12+1 # Number of stress periods
    nstp  = [1]
    for y in range(STRT_YEAR,END_YEAR):
        for m in range(1,13):
            nstp.append(calendar.monthrange(y,m)[1])
    nstp = np.array(nstp)
    steady = np.zeros((nper),dtype=bool)
    steady[0] = True
    
    dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, delr=delr, delc=delc,
                                   top=ztop, botm=botm, perlen=nstp, nstp=nstp, steady=steady, start_datetime='12/31/1983')
        
    # Model Boundaries & initial conditions
    # Active areas
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    for i in range(0,3): ibound[i,:,:] = ibound[i,:,:]*ACTIVE
    ibound[0,:,:] *= np.array(GEO==1)
    
    # Variables for the BAS package
#    strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)
#    strt = H_INIT.get_data(totim=1.0)
    strt = DEM
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
    
    # Layer properties
    # Add LPF package to the MODFLOW model
    # Hydraulic conductivity
    HK_LYR1 = (GEO==1)*ZoneParams[0,0]
    HK_LYR2 = (GEO==1)*ZoneParams[0,1] + (GEO==2)*ZoneParams[0,1] + (GEO==3)*ZoneParams[0,2] \
                + (GEO==4)*ZoneParams[0,3] + (GEO==5)*ZoneParams[0,4]
    HK_LYR3 = ZoneParams[0,5]
    HK = np.array([HK_LYR1,HK_LYR2,HK_LYR3])
    
    # Vertical anisotropy (H:V) of hydraulic conductivity
    VK_LYR1 = (GEO==1)*ZoneParams[1,0]
    VK_LYR2 = (GEO==1)*ZoneParams[1,1] + (GEO==2)*ZoneParams[1,1] + (GEO==3)*ZoneParams[1,2] \
                + (GEO==4)*ZoneParams[1,3] + (GEO==5)*ZoneParams[1,4]
    VK_LYR3 = ZoneParams[1,5]
    VKA = np.array([VK_LYR1,VK_LYR2,VK_LYR3])
    
    # Specific storage
    SS_LYR1 = (GEO==1)*ZoneParams[2,0]
    SS_LYR2 = (GEO==1)*ZoneParams[2,1] + (GEO==2)*ZoneParams[2,1] + (GEO==3)*ZoneParams[2,2] \
                + (GEO==4)*ZoneParams[2,3] + (GEO==5)*ZoneParams[2,4]
    SS_LYR3 = ZoneParams[2,5]
    SS = np.array([SS_LYR1,SS_LYR2,SS_LYR3])
    
    ## Specific yield
    SY_LYR1 = (GEO==1)*ZoneParams[3,0]
    SY_LYR2 = (GEO==1)*ZoneParams[3,1] + (GEO==2)*ZoneParams[3,1] + (GEO==3)*ZoneParams[3,2] \
                + (GEO==4)*ZoneParams[3,3] + (GEO==5)*ZoneParams[3,4]
    SY_LYR3 = ZoneParams[3,5]
    SY = np.array([SY_LYR1,SY_LYR2,SY_LYR3])
    
    lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, hdry=-1e+20, laytyp=[0,0,0], layvka=[1,1,1], 
                                         laywet=[0,0,0], hk=HK, vka=VKA, ss=SS, sy=SY)
    
    return ncol, nrow, mf, dis, bas, lpf

def addNewWells(New_WEL,LYR,WEL_Dict=0,WEL_mult=1,S_YR=0,E_YR=0,coordType='xy',xll=455000,yur=2175000,cellsize=500):
    # New_WEL is an np array of the following format: X (or C), Y (or R), Start Year, End Year, Flow (m3/d)
    # WEL_PAR is a scalar multiplier to be applied to all wells in the data set New_WEL
    # WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for
    #  each well with the layer, row, column, and pumping rate
    # coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array

    # Initialize dictionary    
    if WEL_Dict == 0:
        WEL_Dict = {}
    
    # Assign start year and end year if defined in input
    if S_YR > 0:
        New_WEL[:,2] = np.ones((New_WEL.shape[0]))*S_YR
        
    if E_YR > 0:
        New_WEL[:,3] = np.ones((New_WEL.shape[0]))*E_YR
    
    # Convert X and Y to Column and Row
    if coordType == 'xy':
        cconvert = lambda x: int(np.floor((x - xll)/cellsize))
        New_WEL[:,0] = np.array([cconvert(xi) for xi in New_WEL[:,0]])
        rconvert = lambda y: int(np.floor((yur - y)/cellsize))
        New_WEL[:,1] = np.array([rconvert(yi) for yi in New_WEL[:,1]])
    
    # Loop through all wells in the dataset to fill dictionary
    for w in range(0,New_WEL.shape[0]):
        
        # Assign flow rate for each well to all stress periods indicated by start and end years
        for per in range(int((New_WEL[w,2]-1984)*12+1),int((New_WEL[w,3]-1984)*12+1)):
            try:
                WEL_Dict[per].append([LYR,New_WEL[w,1],New_WEL[w,0],New_WEL[w,4]*WEL_mult])
            except:
                WEL_Dict[per] = [[LYR,New_WEL[w,1],New_WEL[w,0],New_WEL[w,4]*WEL_mult]]
                
    return WEL_Dict

def addRecharge(LU_arrays,PRECIP,S_YR=0,E_YR=0,RCH_Dict=0,RCH_mult=[1,1,1]):
    # Comments**
    
    # Initialize dictionary    
    if RCH_Dict == 0:
        RCH_Dict = {}
    
    if S_YR == 0:
        for l, LU in enumerate(['URBAN','NATURAL','WATER']):
            try:
                RCH_Dict[0] += PRECIP[0]*LU_arrays[l]*RCH_mult[l]
            except:
                RCH_Dict[0] = PRECIP[0]*LU_arrays[l]*RCH_mult[l]
    else:
        for year in range(int(S_YR),int(E_YR)):
            for month in range(1,13):
                
                per = (year-1984)*12+month
                
                for l, LU in enumerate(['URBAN','NATURAL','WATER']):
                    try:
                        RCH_Dict[per] += PRECIP[per]*LU_arrays[l]*RCH_mult[l]
                    except:
                        RCH_Dict[per] = PRECIP[per]*LU_arrays[l]*RCH_mult[l]

    return RCH_Dict

def outputControl(mf):
    # Output Control and Solver
    # Add OC package to the MODFLOW model
    spd = {}
    data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
    spd[0,30] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference' ]
    for m in range(0,12):
        for d in range(0,calendar.monthrange(1984,m+1)[1]):
            spd[m,d] = data2record.copy()
    for y in range(0,30):
        for m in range(1,12):
            spd[y*12+m,(calendar.monthrange(1984+y,m+1)[1]-1)] = data2record.copy()
#    for p in [6,20,32,41,54,78,90,102,114,128,138,150,162,175,187,198,213,225,235]:
#        spd[p,0] = data2record.copy()
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
    
    # Add PCG package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf,mxiter=20, iter1=20)
    
    return oc, pcg

def measureEnergy(hds,WEL_DICT):
    
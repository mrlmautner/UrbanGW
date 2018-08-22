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

def initializeFM(modelname,xll,yll,xur,yur,cellsize,STRT_YEAR,END_YEAR,ACTIVE,GEO,DEM,H_INIT,PhaseParams,ZoneParams):
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
    nper = (END_YEAR - STRT_YEAR)*12 # Number of stress periods
    nstp  = []
    for y in range(STRT_YEAR,END_YEAR):
        for m in range(1,13):
            nstp.append(calendar.monthrange(y,m)[1])
    nstp[0] = 1
    nstp = np.array(nstp)
    steady = np.zeros((nper+1),dtype=bool)
    steady[0] = True
    
    dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, delr=delr, delc=delc,
                                   top=ztop, botm=botm, perlen = nstp, nstp=nstp, steady=steady)
        
    # Model Boundaries & initial conditions
    # Active areas
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    for i in range(0,3): ibound[i,:,:] = ibound[i,:,:]*ACTIVE
    ibound[0,:,:] *= np.array(GEO==1)
    
    # Variables for the BAS package
    strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)
    strt = H_INIT.get_data(totim=1.0)
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=1)
    
    # Layer properties
    # Add LPF package to the MODFLOW model
    # Hydraulic conductivity
    HK_LYR1 = (GEO==1)*PhaseParams[0,0]
    HK_LYR2 = (GEO==1)*PhaseParams[0,1] + (GEO==2)*PhaseParams[0,1] + (GEO==3)*PhaseParams[0,2] 
                + (GEO==4)*PhaseParams[0,3] + (GEO==5)*PhaseParams[0,4]
    HK_LYR3 = PhaseParams[0,5]
    HK = np.array([HK_LYR1,HK_LYR2,HK_LYR3])
    
    # Vertical anisotropy (H:V) of hydraulic conductivity
    VK_LYR1 = (GEO==1)*PhaseParams[1,0]
    VK_LYR2 = (GEO==1)*PhaseParams[1,1] + (GEO==2)*PhaseParams[1,1] + (GEO==3)*PhaseParams[1,2] 
                + (GEO==4)*PhaseParams[1,3] + (GEO==5)*PhaseParams[1,4]
    VK_LYR3 = PhaseParams[1,5]
    VKA = np.array([VK_LYR1,VK_LYR2,VK_LYR3])
    
    # Specific storage
    SS_LYR1 = (GEO==1)*PhaseParams[2,0]
    SS_LYR2 = (GEO==1)*PhaseParams[2,1] + (GEO==2)*PhaseParams[2,1] + (GEO==3)*PhaseParams[2,2] 
                + (GEO==4)*PhaseParams[2,3] + (GEO==5)*PhaseParams[2,4]
    SS_LYR3 = PhaseParams[2,5]
    SS = np.array([SS_LYR1,SS_LYR2,SS_LYR3])
    
    ## Specific yield
    SY_LYR1 = (GEO==1)*PhaseParams[3,0]
    SY_LYR2 = (GEO==1)*PhaseParams[3,1] + (GEO==2)*PhaseParams[3,1] + (GEO==3)*PhaseParams[3,2] 
                + (GEO==4)*PhaseParams[3,3] + (GEO==5)*PhaseParams[3,4]
    SY_LYR3 = PhaseParams[3,5]
    SY = np.array([SY_LYR1,SY_LYR2,SY_LYR3])
    
    lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, hdry=-1e+20, laytyp=[0,0,0], layvka=[1,1,1], 
                                         laywet=[0,0,0], hk=HK, vka=VKA, ss=SS, sy=SY)
    
    return ncol, nrow, mf, dis, bas, lpf

    #%% Wells
    WEL_Dict = {}
    
    # Initialize dictionary with zero fluxes at layer 1, row 1, column 1
    for period in range(0,360):
        WEL_Dict[period] = [[0,1,1,0]]
    
    # Add pumping wells for all years in the study period (extend datasets backwards)
    WEL_Array = np.loadtxt(r'data\RawFiles\PumpingWells\20180430_Pumping_AllDatasets_InModel_WatershedClip_Abridged.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,6,9,10]) # pumping in m3 per day
    for w in range(0,WEL_Array.shape[0]):
        c = int(np.floor((WEL_Array[w,0] - xll)/cellsize))
        r = int(np.floor((yur - WEL_Array[w,1])/cellsize))
        
        for per in range(int(WEL_Array[w,3]),int(WEL_Array[w,4])):
            # Assign pumping multiplier for each parameter period
            if per < PHASE_PER[0]:
                param = 0
            else:
                if PHASE_PER[0] <= per < PHASE_PER[1]:
                    param = 0
                elif PHASE_PER[1] <= per < PHASE_PER[2]:
                    param = 1
                elif PHASE_PER[2] <= per < PHASE_PER[3]:
                    param = 2
                elif PHASE_PER[3] <= per < PHASE_PER[4]:
                    param = 3
                elif PHASE_PER[4] <= per:
                    param = 4
                
                WEL_Dict[per].append([1,r,c,WEL_Array[w,2]*WEL_PAR[param]])
    
    # Add leak wells for all years in the study period
    for param in range(0,4):
        for n in range(0,URBAN[param].shape[0]):
            for per in range(PHASE_PER[param],PHASE_PER[param+1]):
                WEL_Dict[per].append([1,URBAN[param][n,2],URBAN[param][n,3],LEAK_PAR[param]/URBAN[param].shape[0]])
    
    wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_Dict)
    
    #%% Recharge
    RCH_Dict = {}
    precip = []
    
    for year in range(1984,2014):
        for month in range(0,12):
            
            per = (year-1984)*12+month
            RCH_Array = np.zeros((nrow,ncol))
            
            for lu in range(0,3):
                filename = r'data\Input\RCH\ClayMult\RCH' + str(lu+1) + '_' + str(year) + '_' + '{num:02d}'.format(num=month+1) + '.asc'
                precip = openASC(filename)
                RCH_Array += precip*RCH_PAR[lu]
            
            RCH_Dict[per] = RCH_Array
            
    for lu in range(0,3):
        filename = r'data\Input\RCH\RCH' + str(lu+1) + '_AVG.asc'
        precip = openASC(filename)
        RCH_Array += precip*RCH_PAR[lu]
        
        RCH_Dict[0] = RCH_Array
        
    rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_Dict)
    
    #%% Output Control and Solver
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
    for p in [6,20,32,41,54,78,90,102,114,128,138,150,162,175,187,198,213,225,235]:
        spd[p,0] = data2record.copy()
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
    
    # Add PCG package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf)
    
    #%% Run Model and post processing
    # Write the MODFLOW model input files
    print('Data processed in',str(time.time()-timestart),'seconds')
    print('Writing input file...')
    mf.write_input()
    print('Input file written in',str(time.time()-timestart),'seconds')
    
    # Run the MODFLOW model
    success, buff = mf.run_model()
    
    hds = bf.HeadFile(modelname+'.hds')
    h = hds.get_data(totim=1.0)
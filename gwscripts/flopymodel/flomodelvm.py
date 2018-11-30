# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import calendar

def initializeFM(modelname,xll,yll,xur,yur,cellsize,STRT_YEAR,END_YEAR,ACTIVE1,ACTIVE2,GEO,DEM,IH_SS,ZoneParams,THICK1=50,THICK2=350):
    # modelname to set the file root 
    mf = flopy.modflow.Modflow(modelname, exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')
    
    # Model domain and grid definition
    ztop = DEM # Model top elevation (Digital Elevation Model)
    L1botm = ztop - THICK1 # Layer 1 bottom elevation
    L2botm = L1botm - THICK2 # Layer 2 bottom elevation
    botm = [L1botm,L2botm] # Model bottom elevation
    nlay = 2 # Number of layers
    ncol = int((xur-xll)/cellsize) # Number of rows
    nrow = int((yur-yll)/cellsize) # Number of columns
    delr = cellsize # Row height
    delc = cellsize # Column height
    
    # Time discretization
    nper = (END_YEAR - STRT_YEAR)*12 # Number of stress periods
    nstp = []
    for y in range(STRT_YEAR,END_YEAR):
        for m in range(1,13):
            nstp.append(calendar.monthrange(y,m)[1])
    nstp = np.array(nstp)
    steady = np.zeros((nper),dtype=bool)
    
    dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, delr=delr, delc=delc,
                                   top=ztop, botm=botm, perlen=nstp, nstp=nstp, steady=steady, start_datetime='01/01/1984')
        
    # Model Boundaries & initial conditions
    # Active areas
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[0,:,:] = ibound[0,:,:]*ACTIVE1
    ibound[1,:,:] = ibound[1,:,:]*ACTIVE2
    
    # Variables for the BAS package
#    strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)
#    strt = H_INIT.get_data(totim=1.0)
#    strt = DEM
    strt = np.array([IH_SS]*2)
    
    bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
    
    # Layer properties
    # Add LPF package to the MODFLOW model
    # Hydraulic conductivity
    HK_LYR1 = (GEO==1)*ZoneParams[0,0]
    HK_LYR2 = (GEO==1)*ZoneParams[0,1] + (GEO==2)*ZoneParams[0,1] + (GEO==3)*ZoneParams[0,2] \
                + (GEO==4)*ZoneParams[0,3] + (GEO==5)*ZoneParams[0,4]
    HK = np.array([HK_LYR1,HK_LYR2])
    
    # Vertical anisotropy (H:V) of hydraulic conductivity
    VK_LYR1 = (GEO==1)*ZoneParams[1,0]
    VK_LYR2 = (GEO==1)*ZoneParams[1,1] + (GEO==2)*ZoneParams[1,1] + (GEO==3)*ZoneParams[1,2] \
                + (GEO==4)*ZoneParams[1,3] + (GEO==5)*ZoneParams[1,4]
    VKA = np.array([VK_LYR1,VK_LYR2])
    
    # Specific storage
    SS_LYR1 = (GEO==1)*ZoneParams[2,0]
    SS_LYR2 = (GEO==1)*ZoneParams[2,1] + (GEO==2)*ZoneParams[2,1] + (GEO==3)*ZoneParams[2,2] \
                + (GEO==4)*ZoneParams[2,3] + (GEO==5)*ZoneParams[2,4]
    SS = np.array([SS_LYR1,SS_LYR2])
    
    ## Specific yield
    SY_LYR1 = (GEO==1)*ZoneParams[3,0]
    SY_LYR2 = (GEO==1)*ZoneParams[3,1] + (GEO==2)*ZoneParams[3,1] + (GEO==3)*ZoneParams[3,2] \
                + (GEO==4)*ZoneParams[3,3] + (GEO==5)*ZoneParams[3,4]
    SY = np.array([SY_LYR1,SY_LYR2])
    
    lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, hdry=-1e+20, laytyp=[0,0], layvka=[1,1], 
                                         laywet=[0,0], hk=HK, vka=VKA, ss=SS, sy=SY)
    
    return ncol, nrow, mf, dis, bas, lpf

def addNewWells(New_WEL,LYR,mun,WEL_Dict=0,INFO_Dict=0,WEL_mult=1,S_YR=0,E_YR=0,dateType='yr',coordType='xy',xll=455000,yur=2175000,cellsize=500,wellType=0,munleak=1,F=1):
    # New_WEL is an np array of the following format: X (or C), Y (or R), Start Year, End Year, Flow (m3/d)
    # WEL_PAR is a scalar multiplier to be applied to all wells in the data set New_WEL
    # WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for
    #  each well with the layer, row, column, and pumping rate
    # coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array
    # dateType is a marker that should be either 'yr' or 'per' depending on whether a value of year or stress period is passed in the Start and End time columns
    # wellType is a marker that corresponds to the well source/sink type as follows: 0=pumping, 1=good quality (recharge basin,LID surface,WWTP injection), -1=poor quality(leak,LID deep,WWTP surface)
    # munleak is an array with leak % by municipality if the data set is pumping wells and 1 otherwise
    # F is the percentage of leak that remains afer the fix (ie 1 means there have been no fixes)
    
    # Initialize dictionary    
    if WEL_Dict == 0:
        WEL_Dict = {}
    if INFO_Dict == 0:
        INFO_Dict = {}
    
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
    
    if coordType == 'rc':
        New_WEL[:,0] = np.array([int(xi) for xi in New_WEL[:,0]-1])
        New_WEL[:,1] = np.array([int(yi) for yi in New_WEL[:,1]-1])
        
    if dateType == 'yr':
        New_WEL[:,2] = (New_WEL[:,2]-1984)*12+1
        New_WEL[:,3] = (New_WEL[:,3]-1984)*12+1
    
    # Loop through all wells in the dataset to fill dictionary
    for w in range(0,New_WEL.shape[0]):
        r = New_WEL[w,1]
        c = New_WEL[w,0]
        wellmun = mun[int(r),int(c)]
        LEAK_mult = 1
        
        # Reduce the pumping amount by the amount saved by fixing leaks
        if type(munleak) is not int:
            P = float(munleak[np.where(munleak==wellmun)[0],2]) # the percent of the total pumping that is made up of leaks
            LEAK_mult *= (1 - P*(1 - F))
    
        # Assign flow rate for each well to all stress periods indicated by start and end years
        for per in range(int(New_WEL[w,2]-1),int(New_WEL[w,3]-1)):
            try:
                WEL_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult])
                INFO_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult,wellmun,wellType]) # layer, row, column, volume (m3/d), municipality, well type
            except:
                WEL_Dict[per] = [[LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult]]
                INFO_Dict[per]= [[LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult,wellmun,wellType]]
                
    return WEL_Dict,INFO_Dict

def addRecharge(LU_arrays,PRECIP,S_YR=0,E_YR=0,RCH_Dict=0,RCH_mult=[1,1,1]):
    # LU_arrays: dictionary with 3 eantries, one for each land use type which
    #   contains gridded percent amounts for each land use type
    # PRECIP: dictionary with 361 entries, one for each stress period which
    #   contains gridded precipitation
    # RCH_Dict: existing dictionary holding recharge data or 0 if the dictionary
    #   must be initialized
    
    # Initialize dictionary: if there is no exisiting dictionary, create dictionary with no entries
    if RCH_Dict == 0:
        RCH_Dict = {}
    
    # If the recharge is for the first time step, apply only to the first time step
    if S_YR == 0:
        for l, landuse in enumerate(['URBAN','NATURAL','WATER']):
            # If there is not already an entry for the selected stress period, create a new array
            try:
                RCH_Dict[0] += PRECIP[0]*LU_arrays[landuse]*RCH_mult[l]
            except:
                RCH_Dict[0] = PRECIP[0]*LU_arrays[landuse]*RCH_mult[l]
    # Loop through all stress periods between S_YR and E_YR
    else:
        for year in range(int(S_YR),int(E_YR)):
            for month in range(0,12):
                
                per = (year-1984)*12+month
                
                # Apply recharge amounts for each land use type                
                for l, landuse in enumerate(['URBAN','NATURAL','WATER']):                    
                    # If there is not already an entry for the selected stress period, create a new array
                    try:
                        RCH_Dict[per] += PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
                    except:
                        RCH_Dict[per] = PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]

    return RCH_Dict

def outputControl(mf):
    # Output Control and Solver
    # Add OC package to the MODFLOW model
    spd = {}
    data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
    for y in range(0,30):
        for m in range(1,13):
#            spd[y*12+m-1,(calendar.monthrange(1984+y,m)[1]-1)] = data2record.copy()
            for d in range(0,calendar.monthrange(1984+y,m)[1]):
                spd[y*12+m-1,d] = data2record.copy()
    spd[14,30] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference']
#    for p in [6,20,32,41,54,78,90,102,114,128,138,150,162,175,187,198,213,225,235]:
#        spd[p,0] = data2record.copy()
    oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)
#    oc = flopy.modflow.ModflowOc.load('VM_Test.oc', mf,ext_unit_dict=ext_unit_dict)
    
    # Add PCG package to the MODFLOW model
    pcg = flopy.modflow.ModflowPcg(mf,mxiter=20, iter1=20)
    
    return oc, pcg
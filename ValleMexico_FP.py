# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import flopy.utils.binaryfile as bf
import numpy as np
import time
from gwscripts.dataprocessing import gengriddata as gen
from gwscripts.flopymodel import flomodelvm as mod

def run_scenario_model(scenario,num_WWTP,num_RCHBASIN,fixleak,seed=1):
    '''
    num_WWTP is the number of WWTPs to rehabilitate
    num_RCHBASIN is the number of recharge basins to install
    fixleak is the ratio of fixed leaks to historical leaks, 1 indicates same level as historical leaks and 0 indicates all leaks are fixed
    '''
    
    timestart = time.time()
    
    xll = 455000
    yll = 2107000
    xur = 539000
    yur = 2175000
    cellsize = 500
    
    STRT_YEAR = 1984
    END_YEAR = 2014
    
    print('Processing data...')
    
    # Load datasets
    ACTIVE_LYR1 = gen.openASC('data_output\ACTIVE_VM_LYR1.asc') # Extent of model layer 1
    ACTIVE_LYR2 = gen.openASC('data_output\ACTIVE_VM_LYR2.asc') # Extent of model layer 2
    TH1 = gen.openASC('data_output\THICK1_VM.asc') # Thickness of model layer 1
    TH2 = gen.openASC('data_output\THICK2_VM.asc') # Thickness of model layer 2
    GEO = gen.openASC('data_output\GEO_VM.asc') # Geologic formations in layer 2
    DEM = gen.openASC('data_output\DEM_VM.asc') # Digital elevation model of the basin
    IH = gen.openASC('data_output\IH_1984.asc') # Initial hydraulic head in layer 1 and layer 2
    MUN = gen.openASC('data_output\MUN_VM.asc') # Geographic extent of each municipality
    
    # Assign name and create modflow model object
    modelname = 'model_output\VM_' + scenario
    
    ### Parameter data
    
    # Model phase specific parameters
    # Phase starting stress period
    PHASE_PER = [0, 132, 252, 360]
    phases = len(PHASE_PER)-1
    S_per = PHASE_PER[0:len(PHASE_PER)-1]
    E_per = PHASE_PER[1:len(PHASE_PER)]
    # Phase land use dataset year
    LU_PAR = ['1990', '2000', '2010']
    # Phase well pumping multiplier 
    WEL_PAR = np.array([2.671E+00,2.581E+00,2.558E+00])
    # Phase distribution system leak multiplier
    LEAK_PAR = np.array([1,1,1])
    # Percentage of groundwater pumping as ratio of total water use
    GW_to_WU = [0.7131,0.644,0.574966] 
    # Phase LID increase multiplier
    LID_PAR = [1, 1, 1]
    
    # Geologic zone specific parameters translated into an array to input into the model
    # COL1 : Zone hydraulic conductivity vertical anisotropy
    VK_PAR = [100,100,10,1,1]
    # COL2 : Zone specific storage
    SS_PAR = [6.562E-02,1.073E-03,3.176E-02,1.214E-05,9.483E-06]
    # COL3 : Zone horizontal hydraulic conductivity
    HK_PAR = [4.320E-04,2.331E+02,6.750E-02,2.518E-01,8.640E-02] # m/d
    # COL4 : Zone specific yield
    SY_PAR = [0.06,0.15,0.15,0.30,0.01]
    ZoneParams = np.array([HK_PAR,VK_PAR,SS_PAR,SY_PAR])
    
    municipalities = np.unique(MUN)[1:]
    
    RCH_PAR = [1.00E-02, 1.949E-01, 50.00E-02] # Recharge multiplier for urban, natural, and water cover
    
    # Initialize the modflow model with the boundary conditions input above
    ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                     STRT_YEAR,END_YEAR,ACTIVE_LYR1,ACTIVE_LYR2,GEO,DEM,IH,
                                                     ZoneParams,THICK1=TH1,THICK2=TH2)
    
    print('Basic, Discretization, and Layer packages generated in',str(time.time()-timestart),'seconds')
    
    ### Land Use Type
    # Fill a land use dictionary with the ARRAYs that represent the % of each land use cover in each cell
    # and the LISTs that contain all the cells and percentages of each land use type
    LU = {}
    for i, LUset in enumerate(LU_PAR):
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
    
    ### Recharge
    # Create recharge dictionary for MODFLOW RCH package based on land use multipliers and interpolated precipitation rasters
    newtime = time.time()
        
    RCH_DICT = {}
    Precip_Dict = {}
        
    for year in range(int(STRT_YEAR),int(END_YEAR)):
        for month in range(1,13):
            per = (year-STRT_YEAR)*12+month-1
        
            filename = r'data_output\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
            Precip_Dict[per] = gen.openASC(filename)
    
    for i, LUset in enumerate(LU_PAR):
        RCH_DICT = mod.addRecharge(LU_arrays=LU[LUset]['ARRAY'],PRECIP=Precip_Dict,start=S_per[i]+1,end=E_per[i]+1,RCH_Dict=RCH_DICT,RCH_mult=RCH_PAR)
    
    # Create MODFLOW RCH package
    rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)
    
    print('RCH_Dict generated in',str(time.time()-newtime),'seconds')
    
    ### Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
    newtime = time.time()
    LEAK_MUN = np.loadtxt('data_output\leak\LEAK_TOT_MUN.csv',delimiter=',',skiprows=1) # Total recharge percent per municipality: equal to percent of total water use (1997 values) x percent leak (~35%) x recharge percent (15%)
    WWTPs = np.loadtxt(r'data_output\scenarios\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6,11])
    
    # Get groundwater percentage of total based on time period
    gval = np.zeros(360)
    gval[:PHASE_PER[1]] = GW_to_WU[0]
    gval[PHASE_PER[1]:PHASE_PER[2]] = GW_to_WU[1]
    gval[PHASE_PER[2]:] = GW_to_WU[2]
        
    WEL_DICT = {}
    WEL_INFO = {}
    
    # Add supply wells, includes the fixleak and LEAK_MUN terms to reduce pumping when leaks are fixed
    PUMP_PARAM = [WEL_PAR[0],WEL_PAR[1],WEL_PAR[2]]
    
    # Import CONAGUA and SACM pumping datasets
    PUMP_array = np.loadtxt(r'data_output\wells\PUMP_C.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
    WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
    PUMP_array = np.loadtxt(r'data_output\wells\PUMP_S.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
    WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
    
    # Generate monthly pumping datasets for REPDA data in single pumping value format
    for i in range(phases):
        PUMP_array = np.loadtxt(r'data_output\wells\PUMP_RC_Q.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
        
        WEL_DICT, WEL_INFO = mod.addNewWells(New_WEL=PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=PUMP_PARAM[i],start=S_per[i]+1,end=E_per[i]+1,mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
    
    # Sum total pumping
    total_mthly_pumping = np.zeros(PHASE_PER[phases])
    for i in range(total_mthly_pumping.shape[0]):
        total_mthly_pumping[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
    
    # Include only municipalities with urban land cover
    mun = municipalities.copy()
    mun = np.delete(mun,29)
    mun = np.delete(mun,21)
    
    LEAK_arrays = {}
    for i, leakset in enumerate(LU_PAR):
        LEAK_arrays[leakset] = np.zeros((LU[leakset]['LIST']['URBAN'].shape[0]*(PHASE_PER[i+1]-PHASE_PER[i]),5))
        j = 0
        
        # Add leak wells for each period
        for p in range(PHASE_PER[i],PHASE_PER[i+1]):
            for n,m in enumerate(mun):
                # Calculate the total urban area in normalized cell area (1 cell = 1) per municipality
                tempLeak = LU[leakset]['LIST']['URBAN'][(LU[leakset]['LIST']['URBAN'][:,4]==m),:4]
                u_area = tempLeak[:,2].sum()
                
                # Use total pumping for each stress period to determine leak quantities
                LperArea = (-1*total_mthly_pumping[p]/GW_to_WU[i])*fixleak*LEAK_MUN[n,i+1]/u_area # Divide total pumping by GW ratio to get total usage, multiply by leak status percent, multiply by % recharge by municipality (FIX)
                tempLeak[:,2] *= LperArea
                tempLeak[tempLeak[:,3]==1,2] *= 0.1 # apply 90% returns to sewer under clay layer (Geologic formation 1)
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),0] = tempLeak[:,0] # Get rows of all cells of urban land use type from list
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),1] = tempLeak[:,1] # Get columns of all cells of urban land use type from list
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),2] = p # Set the period to the current stress period for all urban cells
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),3] = p+1 # Set the end of the period to the next stress period
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),4] = tempLeak[:,2] # Set the multiplier to the percentage of urban land use type stored in list
                j += tempLeak.shape[0] # Set the new index to the previous index plus the number of cells added
            
        LEAK_arrays[leakset] = LEAK_arrays[leakset][(LEAK_arrays[leakset][:,4]>5),:] # Only include cells that contribute at least 5 m3/day
            
        WEL_DICT, WEL_INFO = mod.addNewWells(LEAK_arrays[leakset],LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=LEAK_PAR[i],coordType='rc',mun=MUN,wellType=1)
    
    total_mthly_leak = np.zeros(PHASE_PER[3])
    for i in range(total_mthly_leak.shape[0]):
        total_mthly_leak[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
    total_mthly_leak = (total_mthly_leak - total_mthly_pumping)
    
    # Add WWTP: compute the difference between the installed WWTP minus the actual treatment quantity (m3/s)
    if num_WWTP>0:
        WWTPs[:,3] = WWTPs[:,2] - WWTPs[:,3] # Colmn 2 is installed treatment capacity and column 3 is actual treatment quantity
        WWTPs = np.insert(WWTPs, 2, PHASE_PER[0]+1, axis=1) # Insert starting period
        WWTPs[:,3] = np.ones(WWTPs.shape[0])*PHASE_PER[phases] # Insert ending period
        WWTPs[WWTPs[:,4]<0.01,4] = 0.01 # For any WWTPs with a difference of less than 0.01 m3/s in capacity, assign an injection of 0.01 m3/s
        
        # Randomly select num_WWTP of WWTPs to improve for recharge
        WWTPs = WWTPs[np.random.choice(WWTPs.shape[0],size=num_WWTP,replace=False),:]
        WEL_DICT, WEL_INFO = mod.addNewWells(WWTPs,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=60*60*24,mun=MUN,wellType=-1) # Mult used to convert to m3/d
        
    ### Recharge Basins
    Basin_Array = np.loadtxt(r'data_output\scenarios\RCH_BASIN.csv', delimiter=',', skiprows=1)
    Basins = np.zeros((num_RCHBASIN,2))
    for b in range(0,num_RCHBASIN):
        randBasin = np.random.randint(0,Basin_Array.shape[0])
        c = int(np.floor((Basin_Array[randBasin,0] - xll)/cellsize))
        r = int(np.floor((yur - Basin_Array[randBasin,1])/cellsize))
        Basins[b,:] = [r,c]
        
        for per in range(0,360):
            # Add injection equal to treatment capacity for each parameter period
            WEL_DICT[per].append([1,r,c,1*86400]) # 1 m3/s recharge basins (35 cfs)
            WEL_INFO[per].append([1,r,c,1*86400,MUN[r,c],1])
    
    wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)
        
    print('WEL_Dict generated in',str(time.time()-newtime),'seconds')
    
    ### Generate output control and solver MODFLOW packages 
    oc, pcg = mod.outputControl(mf)
    
    #    hobs = flopy.modflow.ModflowHob.load('ValleMexicoTR_C.ob_hob', mf)

    ### Run Model and post processing
    # Write the MODFLOW model input files
    print('Data processed in',str(time.time()-timestart),'seconds')
    
    newtime = time.time()
    print('Writing input file...')
    
    mf.write_input()
    print('Input file written in',str(time.time()-newtime),'seconds')
    
    # Run the MODFLOW model
    #success, buff = mf.run_model()
    
    return WWTPs, Basins, total_mthly_leak

# AGU Model Runs
WWTPs, Basins, total_pump_leak = run_scenario_model('Test',0,0,1)

#numscenarios = 4
#scenarioList = ['WWTP','Historical','Leak','Basin'] 
#fixleak = [1,1,0.8,1]
#num_WWTP = [74,0,0,0] # 74
#num_RCHBASIN = [0,0,0,5] # 5
#w = [0]*numscenarios
#b = [0]*numscenarios
#l = [0]*numscenarios
#
#for i in range(numscenarios):
#    w[i],b[i],l[i] = run_scenario_model(scenarioList[i],num_WWTP[i],num_RCHBASIN[i],fixleak[i])
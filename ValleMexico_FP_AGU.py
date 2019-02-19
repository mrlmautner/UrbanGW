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
from gwscripts.optimization import measureobjectives as mo
import pickle
from platypus import Problem, Integer, Real, NSGAII

def runscenariomodel(scenario,num_WWTP,num_RCHBASIN,fixleak,seed=1):
    np.random.seed(seed)
#    timestart = time.time()
#    print('Processing data...')
    cost = 0
    fixleak = fixleak/100 # convert from integer to decimal
    
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
    
    # Assign name and create modflow model object
    modelname = 'model_output\VM_'+scenario
    
    #%% Process parameter data
    # COL1: Phase period starting stress period
    PHASE_PER = [0, 132, 252, 360]
    # COL2 : Phase year
    LU_PAR = [1990, 2000, 2010]
    # COL3 : Phase well pumping multiplier 
    WEL_PAR = np.array([2.671E+00,2.581E+00,2.558E+00])
    #COL3 : Phase distribution system leak multiplier
    LEAK_PAR = np.array([1,1,1])
    GW_to_WU = [0.7131,0.644,0.574966] # Percentage of groundwater pumping of total water use
    # COL4 : Phase LID increase multiplier
    LID_PAR = [1, 1, 1]
    # COL5 :  
    PhaseParams = np.array([PHASE_PER,LU_PAR,WEL_PAR,LEAK_PAR,LID_PAR])
    
    # COL1 : Zone hydraulic conductivity vertical anisotropy
    VK_PAR = [100,100,10,1,1]
    # COL2 : Zone specific storage
    SS_PAR = [6.562E-02,1.073E-03,3.176E-02,1.214E-05,9.483E-06]
    # COL3 : Zone hydraulic conductivity
    HK_PAR = [4.320E-04,2.331E+02,6.750E-02,2.518E-01,8.640E-02] # m/d
    # COL4 : Zone specific yield
    SY_PAR = [0.06,0.15,0.15,0.30,0.01]
    ZoneParams = np.array([HK_PAR,VK_PAR,SS_PAR,SY_PAR])
    
    municipalities = np.unique(MUN)[1:]
    #municipalities = dict.fromkeys(municipalities,{'UAREA':{'1990':{},'2000':{},'2010':{}},'OBJECTIVES':np.zeros(3)})
    
    RCH_PAR = [1.000E-02,5.639E-01,5.000E-01] # Recharge multiplier for urban, natural, and water cover
    
    ncol, nrow, mf, dis, bas, lpf = mod.initializeFM(modelname,xll,yll,xur,yur,cellsize,
                                                     STRT_YEAR,END_YEAR,ACTIVE_LYR1,ACTIVE_LYR2,GEO,DEM,IH,
                                                     ZoneParams,THICK1=TH1,THICK2=TH2)
    
#    print('Basic, Discretization, and Layer packages generated in',str(time.time()-timestart),'seconds')
    
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
#    
#    winfofile = 'model_output\objective_data\LU_'+scenario+'.pickle'
#    with open(winfofile, 'wb') as handle:
#        pickle.dump(LU, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    #%% Recharge
#    newtime = time.time()
    
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
    
#    print('RCH_Dict generated in',str(time.time()-newtime),'seconds')
    
    #%% Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
#    newtime = time.time()
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
    start = [1984,1995,2005]
    end = [1995,2005,2014]
    
    PUMP_array = np.loadtxt(r'data_output\wells\PUMP_C.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
    WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,dateType='per',mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
    PUMP_array = np.loadtxt(r'data_output\wells\PUMP_S.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
    WEL_DICT, WEL_INFO = mod.addNewWells(PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,dateType='per',mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
     
    for i, wellset in enumerate(['PUMP_RC_Q','PUMP_RC_Q','PUMP_RC_Q']):
        PUMP_array = np.loadtxt(r'data_output\wells\\' + wellset + '.csv',
                                                  delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
        
        WEL_DICT, WEL_INFO = mod.addNewWells(New_WEL=PUMP_array,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=PUMP_PARAM[i],S_YR=start[i],E_YR=end[i],mun=MUN,munleak=LEAK_MUN,F=fixleak,G=gval)
    
    # Sum total pumping
    total_mthly_pumping = np.zeros(PHASE_PER[3])
    for i in range(total_mthly_pumping.shape[0]):
        total_mthly_pumping[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
    
    # Add leak wells
    mun = municipalities.copy()
    mun = np.delete(mun,29)
    mun = np.delete(mun,21)

    LEAK_arrays = {'1990':{},'2000':{},'2010':{}}
    for i, leakset in enumerate(['1990','2000','2010']):
        LEAK_arrays[leakset] = np.zeros((LU[leakset]['LIST']['URBAN'].shape[0]*(PHASE_PER[i+1]-PHASE_PER[i]),5))
        j = 0
        
        # Add leak wells for each period
        for p in range(PHASE_PER[i],PHASE_PER[i+1]):
            for n,m in enumerate(mun):
                # Calculate the total urban area in normalized cell area (1 cell = 1) per municipality
                tempLeak = LU[leakset]['LIST']['URBAN'][(LU[leakset]['LIST']['URBAN'][:,4]==m),:4]
                u_area = tempLeak[:,2].sum()
                
                LperArea = (-1*total_mthly_pumping[p]/GW_to_WU[i])*fixleak*LEAK_MUN[n,i+1]/u_area # Divide total pumping by 60% to get total usage, multiply by leak status percent, multiply by % recharge by municipality (FIX)
                tempLeak[:,2] *= LperArea
                tempLeak[tempLeak[:,3]==1,2] *= 0.1 # apply 90% returns to sewer under clay layer
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),0] = tempLeak[:,0]
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),1] = tempLeak[:,1]
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),2] = p
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),3] = p+1
                LEAK_arrays[leakset][j:(j+tempLeak.shape[0]),4] = tempLeak[:,2]
                j += tempLeak.shape[0]
            
        LEAK_arrays[leakset] = LEAK_arrays[leakset][(LEAK_arrays[leakset][:,4]>5),:]
            
        WEL_DICT, WEL_INFO = mod.addNewWells(LEAK_arrays[leakset],LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,dateType='per',coordType='rc',mun=MUN,wellType=1)
    
    total_mthly_leak = np.zeros(PHASE_PER[3])
    for i in range(total_mthly_leak.shape[0]):
        total_mthly_leak[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
    total_mthly_leak = (total_mthly_leak - total_mthly_pumping)
    average_leak = total_mthly_leak.sum()/(360*24*60*60) # convert to m3/s
    cost += np.round(2.0030-average_leak,2)*200
        
    # Add WWTP: compute the difference between the installed WWTP minus the actual treatment quantity (m3/s)
    if num_WWTP>0:
        WWTPs[:,3] = WWTPs[:,2] - WWTPs[:,3]
        WWTPs = np.insert(WWTPs, 2, 1984, axis=1)
        WWTPs[:,3] = np.ones(WWTPs.shape[0])*2014
        WWTPs[WWTPs[:,4]<0.01,4] = 0.01
        # Randomly select num_WWTP of WWTPs to improve for recharge
        WWTPs = WWTPs[np.random.choice(WWTPs.shape[0],size=num_WWTP,replace=False),:]
        WEL_DICT, WEL_INFO = mod.addNewWells(WWTPs,LYR=1,WEL_Dict=WEL_DICT,INFO_Dict=WEL_INFO,WEL_mult=60*60*24,mun=MUN,wellType=-1) # Mult used to convert to m3/d
        
        for i in range(WWTPs.shape[0]):
            if WWTPs[i,5] == 3:
                WWTPs[i,5] = 30
            elif WWTPs[i,5] == 2:
                if np.round(WWTPs[i,4],3) >= 0.1:
                    WWTPs[i,5] = 20
                else:
                    WWTPs[i,5] = 4
            else:
                WWTPs[i,5] = 2
                
        cost += WWTPs[:,5].sum()
        
    #%% Recharge Basins
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
        
    cost += num_RCHBASIN*20
    
    wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)
    
#    print('WEL_Dict generated in',str(time.time()-newtime),'seconds')
#    newtime = time.time()
#    
#    winfofile = 'model_output\objective_data\WEL_INFO_'+scenario+'.pickle'
#    with open(winfofile, 'wb') as handle:
#        pickle.dump(WEL_INFO, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    
#    print('WEL_Dict saved in',str(time.time()-newtime),'seconds')
    

    
    #%%
    oc, pcg = mod.outputControl(mf)
    
#    hobs = flopy.modflow.ModflowHob.load('ValleMexicoTR_C.ob_hob', mf)
    
#    print('Data processed in',str(time.time()-timestart),'seconds')
    
#    newtime = time.time()
#    print('Writing input file...')
    
    # Run Model and post processing
    # Write the MODFLOW model input files 
    mf.write_input()
#    print('Input file written in',str(time.time()-newtime),'seconds')
    
    # Run the MODFLOW model
    success, buff = mf.run_model()
    
    return WWTPs, Basins, total_mthly_leak, cost, WEL_INFO, LU

## AGU Model Runs
#WWTPs, Basins, total_pump_leak, cost, WEL_INFO, LU = runscenariomodel('Test',74,5,0.8)
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
#    w[i],b[i],l[i] = runscenariomodel(scenarioList[i],num_WWTP[i],num_RCHBASIN[i],fixleak[i])


#%% 289 Project
def objective_function(x):
    num_WWTP = x[0]
    num_Basin = x[1]
    fix_leak = x[2]
    
    ACTIVE_LYR1 = gen.openASC('data_output\ACTIVE_VM_LYR1.asc')
    TH1 = gen.openASC('data_output\THICK1_VM.asc')
    DEM = gen.openASC('data_output\DEM_VM.asc')
    
    WWTPs, Basins, total_pump_leak, cost, WEL_INFO, LU = runscenariomodel('Opt',num_WWTP,num_Basin,fix_leak)
    heads = bf.HeadFile('model_output\VM_Opt.hds')
    
    energy = mo.measureEnergy(heads,WEL_INFO,DEM)
    subs_array = mo.measureSubidence(heads,DEM,ACTIVE_LYR1,TH1)
    mound_array = mo.measureMound(heads,DEM,ACTIVE_LYR1,LU,[132,252])
    subs = subs_array[0]/subs_array[1]
    mound = mound_array[0]
    
    return [energy,subs,mound,cost]
#
#nfe = 0
#problem = Problem(3, 4)  # define 3 inputs and 4 objectives (and no constraints)
#problem.directions[:] = Problem.MINIMIZE
#wwtp_int = Integer(0,74)
#basin_int = Integer(0, 10)
#leak_int = Integer(0, 100)
#problem.types[:] = [wwtp_int, basin_int, leak_int]
#problem.function = my_scenarios
#algorithm = NSGAII(problem)
#algorithm.run(200)
#
##%%
### grab the variables (note: we are just taking the ones in the location result[0])
#first_variable = algorithm.result[0].variables[0]
#second_variable = algorithm.result[0].variables[1]
#third_variable = algorithm.result[0].variables[2]
#
#print(wwtp_int.decode(first_variable))
#print(basin_int.decode(second_variable))
#print(leak_int.decode(third_variable))
#
##%%
#print(algorithm.result.variables)
#
##%%
#results_list = []
#variable_list = []
#int_list = [wwtp_int, basin_int, leak_int]
#
#for r, results in enumerate(algorithm.result):
#    results_list.append(results.objectives[:])
#    var_list = []
#    for v,var in enumerate(results.variables):
#        var_list.append(int_list[v].decode(var))
#    variable_list.append(var_list)
#
#print(results_list)
#print(variable_list)
#pickle.dump(results_list, open(r'model_output\opt\objectives_200nfe.pkl', "wb" ))
#pickle.dump(variable_list, open(r'model_output\opt\dvariables_200nfe.pkl', "wb" ))

#%%
with open(r'model_output\opt\objectives_200nfe.pkl', 'rb') as handle:
    results_list = pickle.load(handle)
with open(r'model_output\opt\dvariables_200nfe.pkl', 'rb') as handle:
    variable_list = pickle.load(handle)
#%%
def nondom_sort(solutions):
    objs = solutions.copy()
    num_solutions = len(objs)

    # use a boolean index to keep track of nondominated solns
    keep = np.zeros(num_solutions, dtype = bool)

    for i in range(num_solutions):
        for j in range(num_solutions):
            a = objs[i]
            b = objs[j]
            if np.all(a <= b) & np.any(a < b):
                keep[i] = True
                keep[j] = False

    return keep

nondom_VM = nondom_sort(results_list)
npresults = np.array(results_list)
nondom_results = npresults[nondom_VM]

#%%
#import matplotlib.pyplot as plt
#
### parallel axis
#plt.figure()
#for ppoint in nondom_results:
#    ppoint = (ppoint - nondom_results.min(axis=0)) / (nondom_results.max(axis=0) - nondom_results.min(axis=0))
#    plt.plot(range(4), ppoint, 'steelblue')
#
#plt.gca().set_xticks(range(4))
#plt.gca().set_xticklabels(['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])
#plt.show()
#plt.savefig('parallelaxis_200.png')
#
##%%
## NonDom Model Runs
#npresults = np.array(variable_list)
#nondom_vars = npresults[nondom_VM]
#numscenarios = nondom_vars.shape[0]
#w = [0]*numscenarios
#b = [0]*numscenarios
#l = [0]*numscenarios
#
#for i in range(numscenarios):
#    w[i],b[i],l[i],cost, WEL_INFO, LU = runscenariomodel('Opt',nondom_vars[i,0],nondom_vars[i,1],nondom_vars[i,2])
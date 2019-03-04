# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import time
from gwscripts.dataprocessing import gengriddata as gen
import pickle
import calendar

class model():

    # Initializer / Instance attributes
    def __init__(self, scenario, xll, yll, xur, yur, cellsize, strt_yr, end_yr, ACTIVE_LYR1, ACTIVE_LYR2, TH1, TH2, GEO, DEM, IH, MUN, exe_file = r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe'):
        self.name = scenario # Assign name
        self.xll = xll # X coordinate of the lower left corner
        self.yll = yll # Y coordinate of the lower left corner
        self.xur = xur # X coordinate of the upper right corner
        self.yur = yur # Y coordinate of the upper right corner
        self.cellsize = cellsize # Grid size
        self.ncol = int((self.xur - self.xll) / self.cellsize) # Number of rows
        self.nrow = int((self.yur - self.yll) / self.cellsize) # Number of columns
        self.strt_yr = strt_yr
        self.end_yr = end_yr
        self.actv1 = np.loadtxt(ACTIVE_LYR1,skiprows=6) # Extent of model layer 1
        self.actv2 = np.loadtxt(ACTIVE_LYR2,skiprows=6) # Extent of model layer 2
        self.th1 = np.loadtxt(TH1,skiprows=6) # Thickness of model layer 1
        self.th2 = np.loadtxt(TH2,skiprows=6) # Thickness of model layer 2
        self.geo = np.loadtxt(GEO,skiprows=6) # Geologic formations in layer 2
        self.dem = np.loadtxt(DEM,skiprows=6) # Digital elevation model of the basin (model top)
        self.ih = np.loadtxt(IH,skiprows=6) # Initial hydraulic head in layer 1 and layer 2
        self.mun = np.loadtxt(MUN,skiprows=6) # Geographic extent of each municipality
        self.nlay = 2 # This model only accepts 2 layers
        self.exe = exe_file
    
    def initializeFM(self, ZoneParams):
        # modelname to set the file root 
        mf = flopy.modflow.Modflow('model_output\VM_' + self.name, exe_name=self.exe)
        
        # Model domain and grid definition
        L1botm = self.dem - self.th1 # Layer 1 bottom elevation
        L2botm = L1botm - self.th2 # Layer 2 bottom elevation
        botm = [L1botm,L2botm] # Model bottom elevation
        
        # Time discretization
        nper = (self.end_yr - self.strt_yr)*12 # Number of stress periods
        nstp = []
        for y in range(self.strt_yr,self.end_yr):
            for m in range(1,13):
                nstp.append(calendar.monthrange(y,m)[1])
        nstp = np.array(nstp)
        steady = np.zeros((nper),dtype=bool)
        
        dis = flopy.modflow.ModflowDis(mf, nlay=self.nlay, nrow=self.nrow, ncol=self.ncol, nper=nper, delr=self.cellsize, delc=self.cellsize, top=self.dem, botm=botm, perlen=nstp, nstp=nstp, steady=steady, start_datetime='01/01/1984')
            
        # Model Boundaries & initial conditions
        # Active areas
        ibound = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.int32)
        ibound[0,:,:] = ibound[0,:,:]*self.actv1
        ibound[1,:,:] = ibound[1,:,:]*self.actv2
        
        # Variables for the BAS package
        strt = np.array([self.ih]*2)
        
        bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
        
        # Layer properties
        # Add LPF package to the MODFLOW model
        # Hydraulic conductivity
        HK_LYR1 = (self.geo==1) * ZoneParams[0,0]
        HK_LYR2 = (self.geo==1) * ZoneParams[0,1] + (self.geo==2) * ZoneParams[0,1] + (self.geo==3) * ZoneParams[0,2] + (self.geo==4) * ZoneParams[0,3] + (self.geo==5) * ZoneParams[0,4]
        HK = np.array([HK_LYR1, HK_LYR2])
        
        # Vertical anisotropy (H:V) of hydraulic conductivity
        VK_LYR1 = (self.geo==1) * ZoneParams[1,0]
        VK_LYR2 = (self.geo==1) * ZoneParams[1,1] + (self.geo==2) * ZoneParams[1,1] + (self.geo==3) * ZoneParams[1,2] + (self.geo==4) * ZoneParams[1,3] + (self.geo==5) * ZoneParams[1,4]
        VKA = np.array([VK_LYR1, VK_LYR2])
        
        # Specific storage
        SS_LYR1 = (self.geo==1) * ZoneParams[2,0]
        SS_LYR2 = (self.geo==1) * ZoneParams[2,1] + (self.geo==2) * ZoneParams[2,1] + (self.geo==3) * ZoneParams[2,2] + (self.geo==4) * ZoneParams[2,3] + (self.geo==5) * ZoneParams[2,4]
        SS = np.array([SS_LYR1, SS_LYR2])
        
        ## Specific yield
        SY_LYR1 = (self.geo==1) * ZoneParams[3,0]
        SY_LYR2 = (self.geo==1) * ZoneParams[3,1] + (self.geo==2) * ZoneParams[3,1] + (self.geo==3) * ZoneParams[3,2] + (self.geo==4) * ZoneParams[3,3] + (self.geo==5) * ZoneParams[3,4]
        SY = np.array([SY_LYR1, SY_LYR2])
        
        lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, hdry=-1e+20, laytyp=[0,0], layvka=[1,1], 
                                             laywet=[0,0], hk=HK, vka=VKA, ss=SS, sy=SY)
        
        return mf, dis, bas, lpf
    
    def addNewWells(self, New_WEL, LYR, WEL_Dict=0, INFO_Dict=0, WEL_mult=1, start=0, end=0, dateType='per', coordType='xy', wellType=0, munleak=1, F=1, G=1):
        '''
        New_WEL is an np array of the following format: X (or C), Y (or R), Start Year, End Year, Flow (m3/d)
        WEL_PAR is a scalar multiplier to be applied to all wells in the data set New_WEL
        WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for each well with the layer, row, column, and pumping rate
        coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array
        dateType is a marker that should be either 'yr' or 'per' depending on whether a value of year or stress period is passed in the Start and End time columns
        wellType is a marker that corresponds to the well source/sink type as follows: 0=pumping, 1=good quality (recharge basin,LID surface,WWTP injection), -1=poor quality(leak,LID deep,WWTP surface)
        munleak is an array with leak % by municipality if the data set is pumping wells and 1 otherwise
        F is the percentage of leak that has been fixed (ie 1 means all leaks are fixed)
        '''
        
        # Initialize dictionary    
        if WEL_Dict == 0:
            WEL_Dict = {}
        if INFO_Dict == 0:
            INFO_Dict = {}
        
        # Assign start year and end year if defined in input
        if start > 0:
            New_WEL[:,2] = np.ones((New_WEL.shape[0]))*start
            
        if end > 0:
            New_WEL[:,3] = np.ones((New_WEL.shape[0]))*end
        
        # Convert X and Y to Column and Row
        if coordType == 'xy':
            cconvert = lambda x: int(np.floor((x - self.xll) / self.cellsize))
            New_WEL[:,0] = np.array([cconvert(xi) for xi in New_WEL[:,0]])
            rconvert = lambda y: int(np.floor((self.yur - y) / self.cellsize))
            New_WEL[:,1] = np.array([rconvert(yi) for yi in New_WEL[:,1]])
        
        if coordType == 'rc':
            New_WEL[:,0] = np.array([int(xi) for xi in New_WEL[:,0] - 1])
            New_WEL[:,1] = np.array([int(yi) for yi in New_WEL[:,1] - 1])
        
        # Convert data in year format to stress period format (months)
        if dateType == 'yr':
            New_WEL[:,2] = (New_WEL[:,2] - self.strt_yr) * 12 + 1
            New_WEL[:,3] = (New_WEL[:,3] - self.strt_yr) * 12 + 1
        
        # Loop through all wells in the dataset to fill dictionary
        for w in range(0,New_WEL.shape[0]):
            r = New_WEL[w,1]
            c = New_WEL[w,0]
            wellmun = self.mun[int(r),int(c)]
                    
            # Reduce the pumping amount by the amount saved by fixing leaks
            if type(munleak) is not int:
                P = float(munleak[np.where(munleak==wellmun)[0],4]) # the percent of the total usage that is made up of leaks
                                
            # Assign flow rate for each well to all stress periods indicated by start and end years
            for per in range(int(New_WEL[w,2] - 1),int(New_WEL[w,3] - 1)):
                # Determine whether the dataset is a pumping set or any other type
                if type(munleak) is not int:
                    LEAK_mult = 1 - (1 / G[per]) * P * F # Apply a multiplier that subtracts the leak averted from the total water use from the groundwater pumping
                else:
                    LEAK_mult = 1
                
                try:
                    WEL_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult])
                    INFO_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult,wellmun,wellType]) # layer, row, column, volume (m3/d), municipality, well type
                except:
                    WEL_Dict[per] = [[LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult]]
                    INFO_Dict[per]= [[LYR,r,c,New_WEL[w,4]*WEL_mult*LEAK_mult,wellmun,wellType]]
                    
        return WEL_Dict,INFO_Dict
    
    def addRecharge(self,LU_arrays, PRECIP, start=0, end=0, RCH_Dict=0, RCH_mult=[1,1,1], dateType='per'):
        '''
        Outputs a dictionary of recharge arrays based on land use multiplier, land use cover, and precipitation input
        LU_arrays: dictionary with 3 eantries, one for each land use type which contains gridded percent amounts for each land use type
        PRECIP: dictionary with 361 entries, one for each stress period which contains gridded precipitation
        RCH_Dict: existing dictionary holding recharge data or 0 if the dictionary must be initialized dateType: the date format for the start and end variables
        '''
        
        # Initialize dictionary: if there is no exisiting dictionary, create dictionary with no entries
        if RCH_Dict == 0:
            RCH_Dict = {}
        
        # If the recharge is for the first time step, apply only to the first time step
        if start == 0:
            for l, landuse in enumerate(['URBAN','NATURAL','WATER']):
                # If there is not already an entry for the selected stress period, create a new array
                try:
                    RCH_Dict[0] += PRECIP[0] * LU_arrays[landuse] * RCH_mult[l]
                except:
                    RCH_Dict[0] = PRECIP[0] * LU_arrays[landuse] * RCH_mult[l]
        
        # Convert data in year format to stress period format (months)
        if dateType == 'yr':
            start = (start - self.strt_yr) * 12 + 1
            end = (end - self.strt_yr) * 12 + 1
        
        # Loop through all stress periods between S_YR and E_YR
        else:
            for per in range(int(start - 1), int(end - 1)):
                
                # Apply recharge amounts for each land use type                
                for l, landuse in enumerate(['URBAN','NATURAL','WATER']):                    
                    # If there is not already an entry for the selected stress period, create a new array
                    try:
                        RCH_Dict[per] += PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
                    except:
                        RCH_Dict[per] = PRECIP[per]*LU_arrays[landuse]*RCH_mult[l]
    
        return RCH_Dict
    
    def outputControl(self,mf):
        ''' 
        Generate Output Control and Solver packages
        Add OC package to the MODFLOW model
        '''
        spd = {}
        data2record = ['save head', 'save drawdown', 'save budget', 'print budget']
        for y in range(0,30):
            for m in range(1,13):
                for d in range(0,calendar.monthrange(self.strt_yr + y, m)[1]):
                    spd[y * 12 + m - 1, d] = data2record.copy()
        spd[14,30] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference']
        oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

        # Add PCG package to the MODFLOW model
        pcg = flopy.modflow.ModflowPcg(mf,mxiter=20, iter1=20)
        
        return oc, pcg
    
    def run_scenario_model(self,num_WWTP,num_RCHBASIN,fixleak,seed=1):
        '''
        num_WWTP is the number of wastewater treatment plants to rehabilitate for wastewater injection into the aquifer
        num_RCHBASIN is the number of infiltration basins that will recharge the aquifer using imported water
        fixleak is the percent of fixed leaks to historical leaks, 0 indicates the same level as historical leaks and 100 indicates all leaks are fixed
        '''
        
        timestart = time.time()
        print('Processing data...')
        ### Parameter data
        
        # Model phase specific parameters
        # Phase starting stress period
        PHASE_PER = [0, 132, 252, 360]
        phases = len(PHASE_PER) - 1
        S_per = PHASE_PER[0:len(PHASE_PER) - 1]
        E_per = PHASE_PER[1:len(PHASE_PER)]
        # Phase land use dataset year
        LU_PAR = ['1990', '2000', '2010']
        # Phase well pumping multiplier 
        WEL_PAR = np.array([2.671E+00,2.581E+00,2.558E+00])
        # Phase distribution system leak multiplier
        LEAK_PAR = np.array([1,1,1])
        fixleak = fixleak/100 # convert from integer to decimal
        # Percentage of groundwater pumping as ratio of total water use
        GW_to_WU = [0.7131,0.644,0.574966] 
        # Phase LID increase multiplier
        LID_PAR = [1, 1, 1]
        # Initial cost
        cost = 0
        
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
        
        municipalities = np.unique(self.mun)[1:]
        
        RCH_PAR = [1.00E-02, 1.949E-01, 50.00E-02] # Recharge multiplier for urban, natural, and water cover
        
        # Initialize the modflow model with the boundary conditions input above
        mf, dis, bas, lpf = self.initializeFM(ZoneParams)
        
        print('Basic, Discretization, and Layer packages generated in', str(time.time() - timestart), 'seconds')
        
        '''
        Land Use Type
        Fill a land use dictionary with the ARRAYs that represent the % of each land use cover in each cell and the LISTs that contain all the cells and percentages of each land use type
        '''
        LU = {}
        for i, LUset in enumerate(LU_PAR):
            LU[LUset] = {'ARRAY':{},'LIST':{}}
            
            for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
                filename = r'data_output\landuse\LU-' + LUset + '-' + LUtype + '.asc'
                perarea =  np.loadtxt(filename,skiprows=6)
                LU[LUset]['ARRAY'][LUtype] = perarea
                
                LU[LUset]['LIST'][LUtype] = np.zeros((perarea.shape[0]*perarea.shape[1],5))
                
                l = 0
                for row in range(0,perarea.shape[0]):
                    for col in range(0,perarea.shape[1]):
                        if perarea[row,col] > 0.001:
                            LU[LUset]['LIST'][LUtype][l,2] = perarea[row,col]
                        LU[LUset]['LIST'][LUtype][l,0] = col
                        LU[LUset]['LIST'][LUtype][l,1] = row
                        LU[LUset]['LIST'][LUtype][l,3] = 1 - self.actv1[row,col] # 0 if clay layer, 1 if no clay layer
                        LU[LUset]['LIST'][LUtype][l,4] = self.mun[row,col]
                        l += 1
                LU[LUset]['LIST'][LUtype] = LU[LUset]['LIST'][LUtype][LU[LUset]['LIST'][LUtype][:,2]>0,:]
    
        # Save land use database for use in mounding objective
        winfofile = 'model_output\objective_data\LU_' + self.name + '.pickle'
        with open(winfofile, 'wb') as handle:
            pickle.dump(LU, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        '''
        Recharge
        Create recharge dictionary for MODFLOW RCH package based on land use multipliers and interpolated precipitation rasters
        '''
        newtime = time.time()
            
        RCH_DICT = {}
        Precip_Dict = {}
            
        for year in range(int(self.strt_yr),int(self.end_yr)):
            for month in range(1,13):
                per = (year - self.strt_yr) * 12 + month - 1
            
                filename = r'data_output\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
                Precip_Dict[per] = np.loadtxt(filename,skiprows=6)
        
        for i, LUset in enumerate(LU_PAR):
            RCH_DICT = self.addRecharge(LU_arrays=LU[LUset]['ARRAY'], PRECIP=Precip_Dict, start=S_per[i] + 1, end=E_per[i] + 1, RCH_Dict=RCH_DICT, RCH_mult=RCH_PAR)
        
        # Create MODFLOW RCH package
        rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)
        
        print('RCH_Dict generated in', str(time.time() - newtime), 'seconds')
        
        '''
        Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
        '''
        newtime = time.time()
        LEAK_MUN = np.loadtxt('data_output\leak\LEAK_TOT_MUN.csv',delimiter=',',skiprows=1) # Total recharge percent per municipality: equal to percent of total water use (1997 values) x percent leak (~35%) x recharge percent (15%)
        
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
        WEL_DICT, WEL_INFO = self.addNewWells(PUMP_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, munleak=LEAK_MUN, F=fixleak, G=gval)
        PUMP_array = np.loadtxt(r'data_output\wells\PUMP_S.csv',
                                                      delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
        WEL_DICT, WEL_INFO = self.addNewWells(PUMP_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, munleak=LEAK_MUN, F=fixleak, G=gval)
        
        # Generate monthly pumping datasets for REPDA data in single pumping value format
        for i in range(phases):
            PUMP_array = np.loadtxt(r'data_output\wells\PUMP_RC_Q.csv', delimiter=',', skiprows=1, usecols=[1,2,4,5,11]) # pumping in m3 per day
            
            WEL_DICT, WEL_INFO = self.addNewWells(New_WEL=PUMP_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=PUMP_PARAM[i], start=S_per[i] + 1, end=E_per[i] + 1, munleak=LEAK_MUN, F=fixleak, G=gval)
        
        # Sum total pumping
        total_mthly_pumping = np.zeros(PHASE_PER[phases])
        for i in range(total_mthly_pumping.shape[0]):
            total_mthly_pumping[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
        
        # Include only municipalities with urban land cover
        mun = municipalities.copy()
        mun = np.delete(mun,29)
        mun = np.delete(mun,21)
        
        '''
        Leak Repair
        
        Create a well dictionary for all the leak cells. The leak cells will be treated
        as an injection well at each cell in which there is urban land cover. MODFLOW
        distributes injection wells evenly across the area of the cell. The leak
        percentage is based on the leak percentage determined by municipality. Then
        the leak amount determined for each cell is multiplied by the percent of
        urban land cover in that cell. Finally, leaks in cells that are located in the
        lacustrine zone are reduced by 90% assuming that the low hydraulic
        conductivity does not allow for high levels of infiltration and the sewer
        provides a preferential flow path out of the basin
        '''
        LEAK_arrays = {}
        for i, leakset in enumerate(LU_PAR):
            LEAK_arrays[leakset] = np.zeros((LU[leakset]['LIST']['URBAN'].shape[0] * (PHASE_PER[i + 1] - PHASE_PER[i]),5))
            j = 0
            
            # Loop through each period model period
            for p in range(PHASE_PER[i] + 1,PHASE_PER[i + 1] + 1):
                for n,m in enumerate(mun):
                    # Calculate the total urban area in normalized cell area (1 cell = 1) per municipality
                    tempLeak = LU[leakset]['LIST']['URBAN'][(LU[leakset]['LIST']['URBAN'][:,4]==m), :4]
                    u_area = tempLeak[:, 2].sum()
                    
                    # Use total pumping for each stress period to determine leak quantities
                    LperArea = (-1 * total_mthly_pumping[p - 1] / GW_to_WU[i]) * (1 - fixleak) * LEAK_MUN[n, i + 1] / u_area # Divide total pumping by GW ratio to get total usage, multiply by leak status percent, multiply by % recharge by municipality (FIX)
                    tempLeak[:,2] *= LperArea
                    
                    # apply 90% returns to sewer under clay layer (Geologic formation 1)
                    tempLeak[tempLeak[:,3]==1, 2] *= 0.1
                    
                    # Get rows of all cells of urban land use type from list
                    LEAK_arrays[leakset][j:(j + tempLeak.shape[0]), 0] = tempLeak[:, 0]
                    
                    # Get columns of all cells of urban land use type from list
                    LEAK_arrays[leakset][j:(j + tempLeak.shape[0]), 1] = tempLeak[:, 1]
                    
                    # Set the period to the current stress period for all urban cells
                    LEAK_arrays[leakset][j:(j + tempLeak.shape[0]), 2] = p
                    
                    # Set the end of the period to the next stress period
                    LEAK_arrays[leakset][j:(j + tempLeak.shape[0]), 3] = p + 1
                    
                    # Set the multiplier to the percentage of urban land use type stored in list
                    LEAK_arrays[leakset][j:(j + tempLeak.shape[0]), 4] = tempLeak[:,2]
                    
                    # Set the new index to the previous index plus the number of cells added
                    j += tempLeak.shape[0]
                
            LEAK_arrays[leakset] = LEAK_arrays[leakset][(LEAK_arrays[leakset][:, 4] > 5), :] # Only include cells that contribute at least 5 m3/day
                
            WEL_DICT, WEL_INFO = self.addNewWells(LEAK_arrays[leakset], LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=LEAK_PAR[i], coordType='rc', wellType=1)
        
        total_mthly_leak = np.zeros(PHASE_PER[3])
        for i in range(total_mthly_leak.shape[0]):
            total_mthly_leak[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
        total_mthly_leak = (total_mthly_leak - total_mthly_pumping)
        ### Cost attributed to fixing leaks
        average_leak = total_mthly_leak.sum() / (360 * 24 * 60 * 60) # convert to m3/s
        cost += np.round(2.0030 - average_leak, 2)*200
        
        '''
        Wastewater Treatment Reuse
        
        WWTP is imported from a csv that has the following columns:(m3/s)
        The difference between the installed capacity and the actual treatment
        quantity is computed as the injection quantity
        The wastewater treatment plants chosen for rehabilitation are randomly
        selected from the list based on the number of plants indicated
        At the end of the data processing WWTP has the following columns: X, Y, 
        start year, end year, difference between installed and actual capacity
        '''
        WWTPs = np.loadtxt(r'data_output\scenarios\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6,11])
        
        if num_WWTP > 0:
            WWTPs[:, 3] = WWTPs[:,2] - WWTPs[:,3] # Col 2 is installed treatment capacity and Col 3 is actual treatment quantity
            WWTPs = np.insert(WWTPs, 2, PHASE_PER[0] + 1, axis=1) # Insert starting period
            WWTPs[:, 3] = np.ones(WWTPs.shape[0]) * PHASE_PER[phases] # Replace Col 3 with ending period
            WWTPs[WWTPs[:, 4] < 0.01, 4] = 0.01 # For any WWTPs with a difference of
            # less than 0.01 m3/s in capacity, assign an injection of 0.01 m3/s
            
            WWTPs = WWTPs[np.random.choice(WWTPs.shape[0], size=num_WWTP, replace=False), :]
            WEL_DICT, WEL_INFO = self.addNewWells(WWTPs, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=60*60*24, wellType=-1) # Mult used to convert to m3/d
            
            ## Cost attributed to each exapanded WWTPs
            for i in range(WWTPs.shape[0]):
                if WWTPs[i, 5] == 3:
                    WWTPs[i, 5] = 30
                elif WWTPs[i, 5] == 2:
                    if np.round(WWTPs[i, 4], 3) >= 0.1:
                        WWTPs[i, 5] = 20
                    else:
                        WWTPs[i, 5] = 4
                else:
                    WWTPs[i, 5] = 2
                    
            cost += WWTPs[:, 5].sum()
        else:
            WWTPs = np.zeros(5)
        
        '''
        Recharge Basins
        '''
        Basin_Array = np.loadtxt(r'data_output\scenarios\RCH_BASIN.csv', delimiter=',', skiprows=1)
        Basins = np.zeros((num_RCHBASIN, 2))
        for b in range(0, num_RCHBASIN):
            randBasin = np.random.randint(0, Basin_Array.shape[0])
            c = int(np.floor((Basin_Array[randBasin,0] - self.xll) / self.cellsize))
            r = int(np.floor((self.yur - Basin_Array[randBasin, 1]) / self.cellsize))
            Basins[b, :] = [r, c]
            
            for per in range(0, 360):
                # Add injection equal to treatment capacity for each parameter period
                WEL_DICT[per].append([1, r, c, 1 * 86400]) # 1 m3/s recharge basins (35 cfs)
                WEL_INFO[per].append([1, r, c, 1 * 86400, self.mun[r, c], 1])
        
        ## Cost attributed to building recharge basins
        cost += num_RCHBASIN * 20
        
        wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT)
            
        print('WEL_Dict generated in', str(time.time() - newtime), 'seconds')
    
        ## Create pickle file of Well Data to be available for post processing of well energy use objective
        winfofile = 'model_output\objective_data\WEL_INFO_' + self.name + '.pickle'
        with open(winfofile, 'wb') as handle:
            pickle.dump(WEL_INFO, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
        print('WEL_Dict saved in',str(time.time()-newtime),'seconds')
        
        # Generate output control and solver MODFLOW packages 
        oc, pcg = self.outputControl(mf)
        
        ##    hobs = flopy.modflow.ModflowHob.load('ValleMexicoTR_C.ob_hob', mf)
    
        # Run Model and post processing
        ## Write the MODFLOW model input files
        print('Data processed in', str(time.time() - timestart), 'seconds')
        
        newtime = time.time()
        print('Writing input file...')
        
        mf.write_input()
        print('Input file written in', str(time.time() - newtime), 'seconds')
        
        newtime = time.time()
        print('Running MODFLOW model...')
        # Run the MODFLOW model
        success, buff = mf.run_model(silent=True)
        
        print('MODFLOW model completed run in', str(time.time() - newtime), 'seconds')
        
        self.wwtps = WWTPs
        self.basins = Basins
        self.mthlyleak = total_mthly_leak
        self.cost = cost
        self.wells = WEL_INFO
        self.landuse = LU

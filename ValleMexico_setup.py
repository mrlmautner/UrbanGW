# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:23:15 2018

@author: MM
"""

import flopy
import numpy as np
import time
import pickle
import calendar

class model():

    # Initializer / Instance attributes
    def __init__(self, scenario, xll, yll, xur, yur, cellsize, strt_yr, end_yr, ACTIVE, THICKNESS, GEO, DEM, IH, MUN, PAR, exe_file = r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe'):
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
        self.actv = [np.loadtxt(i,skiprows=6) for i in ACTIVE]  # Extent of model layers 1 through n
        self.thck = [np.loadtxt(i,skiprows=6) for i in THICKNESS] # Thickness of model layers 1 through n
        self.geo = [np.loadtxt(i,skiprows=6) for i in GEO] # Geologic formations in layers 1 through n
        self.dem = np.loadtxt(DEM,skiprows=6) # Digital elevation model of the basin (model top)
        self.ih = np.loadtxt(IH,skiprows=6) # Initial hydraulic head in layer 1 and layer 2
        self.mun = np.loadtxt(MUN,skiprows=6) # Geographic extent of each municipality
        self.nlay = 2 # This model only accepts 2 layers
        self.exe = exe_file
        
        # Create adjustable parameter dictionary
        self.params = {}
        with open(PAR) as f:
            pval = f.readlines()
        
        for i in pval:
            try:
                self.params[i[:i.rfind('_')]].append(float(i[i.rfind('   '):]))
            except:
                try:
                    self.params[i[:i.rfind('_')]] = [float(i[i.rfind('   '):])]
                except:
                    pass
    
    def initializeFM(self):
        # modelname to set the file root 
        mf = flopy.modflow.Modflow(r'model_files\modflow\VM_' + self.name, exe_name=self.exe)
        
        # Model domain and grid definition
        botm = np.zeros((self.nlay, self.nrow, self.ncol))
        sumthck = np.zeros((self.nrow, self.ncol))
        for b, thickness in enumerate(self.thck):
            sumthck = np.add(sumthck,thickness)
            botm[b,:,:] = self.dem - sumthck # Current layer bottom elevation
        self.botm = botm
        
        # Time discretization
        nper = (self.end_yr - self.strt_yr)*12 # Number of stress periods
        nstp = []
        for y in range(self.strt_yr,self.end_yr):
            for m in range(1,13):
                nstp.append(calendar.monthrange(y,m)[1])
        nstp = np.array(nstp)
        steady = np.zeros((nper),dtype=bool)
        
        dis = flopy.modflow.ModflowDis(mf, nlay=self.nlay, nrow=self.nrow, ncol=self.ncol, nper=nper, delr=self.cellsize, delc=self.cellsize, top=self.dem, botm=botm, perlen=nstp, nstp=9, tsmult=1.3, steady=steady, start_datetime='01/01/1984')
            
        # Model Boundaries & initial conditions
        # Active areas
        ibound = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.int32)
        for a, active in enumerate(self.actv):
            ibound[a,:,:] = active # Current layer active area
        
        # Variables for the BAS package
        strt = np.array([self.ih]*2)
        
        bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=3)
        
        # Layer properties
        # Add LPF package to the MODFLOW model
        # Create a dictionary of arrays with geologic characteristics: HK = Hydraulic conductivity, VANI = Vertical anisotropy (H:V) of hydraulic conductivity, SS = Specific storage, SY = Specific yield
        geoarrays = {}
        
        # Loop through the layers and formations for each layer to apply the geologic parameters to each array
        for p in ['HK', 'VANI', 'SS', 'SY']:   
            geoarrays[p] = np.zeros((self.nlay,self.nrow,self.ncol))
            
            for l in range(self.nlay):
                for f, fval in enumerate(self.params[p]):
                    geoarrays[p][l,:,:] += (self.geo[l] == f+1) * fval
            
        layvka = [1]*self.nlay # Indicates that VANI represents the ratio of H:V hydraulic conductivity
        
        lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, laytyp=[1,1], layvka=layvka, hk=geoarrays['HK'], vka=geoarrays['VANI'], ss=geoarrays['SS'], sy=geoarrays['SY'], laywet=[1,1])
#        lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, laytyp=[0,0], layvka=layvka, hk=geoarrays['HK'], vka=geoarrays['VANI'], ss=geoarrays['SS'], sy=geoarrays['SY'])
        
        return mf, dis, bas, lpf
    
    def addNewWells(self, New_WEL, LYR, WEL_Dict=0, INFO_Dict=0, WEL_mult=1, start=0, end=0, dateType='per', coordType='xy', pumpwell=False, changepumping=False):
        '''
        New_WEL is an np array of the following format: X (or C), Y (or R), Start Year, End Year, Flow (m3/d)
        WEL_mult is a scalar multiplier to be applied to all wells in the data set New_WEL
        WEL_Dict is a dictionary that contains dictionary for each stress period, each dictionary contains an entry for each well with the layer, row, column, and pumping rate
        coordType is a marker that should be either 'xy' or 'rc' depending on the coordinate definition in the New_WEL array
        dateType is a marker that should be either 'yr' or 'per' depending on whether a value of year or stress period is passed in the Start and End time columns
        pumpwell
        changepumping
        '''
        
        # Initialize dictionary    
        if WEL_Dict == 0:
            WEL_Dict = {}
        if INFO_Dict == 0:
            INFO_Dict = {}
        
        # Assign start period and end period if defined in input
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
                    
            # Reduce the pumping amount by a percentage by municipality
            if changepumping:
                P = float(self.altpump[np.where(self.altpump==wellmun)[0],1]) # the ratio of new pumping to old pumping
            else:
                P = 1
                    
            # Assign flow rate for each well to all stress periods indicated by start and end years
            for per in range(int(New_WEL[w,2] - 1),int(New_WEL[w,3] - 1)):
                if pumpwell:
                    R = self.ratiogn['PER'][per]
                else:
                    R = 1
                
                try:
                    WEL_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*P*R])
                    INFO_Dict[per].append([LYR,r,c,New_WEL[w,4]*WEL_mult*P*R,wellmun]) # layer, row, column, volume (m3/d), municipality, well type
                except:
                    WEL_Dict[per] = [[LYR,r,c,New_WEL[w,4]*WEL_mult*P*R]]
                    INFO_Dict[per]= [[LYR,r,c,New_WEL[w,4]*WEL_mult*P*R,wellmun]]
                    
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
                spd[y * 12 + m - 1, 8] = data2record.copy()
#                spd[y * 12 + m - 1, calendar.monthrange(self.strt_yr + y, m)[1] - 1] = data2record.copy() # If time steps in month is equal to number of days
#                for d in range(0,calendar.monthrange(self.strt_yr + y, m)[1]):
#                    spd[y * 12 + m - 1, d] = data2record.copy()
        spd[26,8] = ['save head', 'save drawdown', 'save budget', 'print budget', 'ddreference']
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
        
        np.random.seed(seed)
        timestart = time.time()
        print('Processing data...')
        
        # Phase starting stress period
        PHASE_PER = [0, 132, 252, 360]
        phases = len(PHASE_PER) - 1
        S_per = PHASE_PER[0:len(PHASE_PER) - 1]
        E_per = PHASE_PER[1:len(PHASE_PER)]
        # Phase land use dataset year
        LU_PAR = ['1990', '2000', '2010']
        
        # Model internal variables
        drains = True
        fixleak = fixleak/100 # convert from integer to decimal
        cost = 0 # Initial cost
        sec2day = 60*60*24 # Second to day conversion
        LID_PAR = [1, 1, 1] # Phase LID increase multiplier
        
        # Water supply data
        hist_water_use = np.loadtxt(r'model_files\optimization_data\decisions\twu.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Initial (original) all other supplies before alternatives, matrix of size municipalities by phases (m3/s)
        total_water_use = hist_water_use*self.params['TWU']*sec2day # Multiply by total water use parameters (m3/d)
        self.twateruse = total_water_use.sum(axis=0) # Total water use for each model phase (m3/d)
        
        i_other = np.loadtxt(r'model_files\optimization_data\decisions\initial_supply.csv', delimiter=',', skiprows=1, usecols=(1,2,3)) # Initial (original) all other supplies before alternatives (m3/s)
        i_other = i_other.sum(axis=0)*sec2day # Initial other water supply for each model phase (m3/d)
        new_other = i_other.copy() # New other water supply for each model phase (m3/d)
        
        # Alternatives changes to supply
        # Import alternative pumping scheme (percentage changes in pumping), total groundwater pumping must be equal to original
        self.altpump = np.loadtxt(r'model_files\optimization_data\decisions\altpump.csv', delimiter=',', skiprows=1)
        
        # Calculate historical quantity of leaks in each municipality
        LEAK_MUN = np.loadtxt(r'data_processed\leak\LEAK_TOT_MUN.csv',delimiter=',',skiprows=1) # Total recharge percent per municipality: equal to percent of total water use (1997 values) x percent leak (~35%) x recharge percent (15%)
        leaks = np.zeros((LEAK_MUN.shape[0],phases+1))
        leaks[:,0] = LEAK_MUN[:,0]
        for i in range(phases):
            leaks[:,i+1] = self.params['LK'][i]*LEAK_MUN[:,1]*total_water_use[:,i] # Total leak per municipality by model phase (m3/d)
        
        new_other += leaks[:,1:].sum(axis=0) * fixleak # Add the leaks averted as an alternative supply
        self.ratiogn = {}
        self.ratiogn['PHASE'] = (self.twateruse - new_other)/(self.twateruse - i_other) # Create a ratio of groundwater use with new other supplies to groundwater use with initial other supplies to apply to all groundwater pumping (dimensionless)
        self.ratiogn['PER'] = np.zeros(PHASE_PER[phases])
        for i, LUset in enumerate(LU_PAR):
            for p in range(PHASE_PER[i],PHASE_PER[i+1]):
                self.ratiogn['PER'][p] = self.ratiogn['PHASE'][i]
        
        # Initialize the modflow model with the boundary conditions input above
        mf, dis, bas, lpf = self.initializeFM()
        
        print('Basic, Discretization, and Layer packages generated in', str(time.time() - timestart), 'seconds')
        
        '''
        Land Use Type
        Fill a land use dictionary with the ARRAYs that represent the % of each land use cover in each cell and the LISTs that contain all the cells and percentages of each land use type
        '''
        LU = {}
        for i, LUset in enumerate(LU_PAR):
            LU[LUset] = {'ARRAY':{},'LIST':{}}
            
            for l, LUtype in enumerate(['URBAN','NATURAL','WATER']):
                filename = r'data_processed\landuse\LU-' + LUset + '-' + LUtype + '.asc'
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
                        LU[LUset]['LIST'][LUtype][l,3] = 1 - self.geo[0][row,col] # 0 if clay layer, 1 if no clay layer
                        LU[LUset]['LIST'][LUtype][l,4] = self.mun[row,col]
                        l += 1
                LU[LUset]['LIST'][LUtype] = LU[LUset]['LIST'][LUtype][LU[LUset]['LIST'][LUtype][:,2]>0,:]
    
        # Save land use database for use in mounding objective
        winfofile = r'model_files\optimization_data\objectives\LU_' + self.name + '.pickle'
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
            
                filename = r'data_processed\recharge\claymult\PrecipCM_' + str(year) + '_' + '{num:02d}'.format(num=month) + '.asc'
                Precip_Dict[per] = np.loadtxt(filename,skiprows=6)
        
        for i, LUset in enumerate(LU_PAR):
            RCH_DICT = self.addRecharge(LU_arrays=LU[LUset]['ARRAY'], PRECIP=Precip_Dict, start=S_per[i]+1, end=E_per[i]+1, RCH_Dict=RCH_DICT, RCH_mult=self.params['RCH'])
        
        # Create MODFLOW RCH package
        rch = flopy.modflow.ModflowRch(mf, nrchop=3,  ipakcb=9, rech=RCH_DICT)
        
        print('RCH_Dict generated in', str(time.time() - newtime), 'seconds')
        
        '''
        Well objects: supply wells, distribution leaks, injection wells, wastewater reuse, recharge basins
        '''
        newtime = time.time()
        
        WEL_DICT = {}
        WEL_INFO = {}
        
        # Add supply wells, includes the ratioGn multiplier to reduce pumping when new supplies are added
        # Import CONAGUA and SACM pumping datasets
        CAEM_array = np.loadtxt(r'data_processed\wells\PUMP_C.csv', delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
        WEL_DICT, WEL_INFO = self.addNewWells(CAEM_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, pumpwell=True)
            
        SACM_array = np.loadtxt(r'data_processed\wells\PUMP_S-ERRORS.csv',delimiter=',', skiprows=1, usecols=[1,2,7,8,11]) # pumping in m3 per day
        WEL_DICT, WEL_INFO = self.addNewWells(SACM_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, pumpwell=True)
        
        REPDA_array = np.loadtxt(r'data_processed\wells\PUMP_RC_Q-Repeats.csv', delimiter=',', skiprows=1, usecols=[1,2,4,5,11,16]) # pumping in m3 per day
        total_urban_repda = REPDA_array[REPDA_array[:,5]==1,4].sum() * -1 # Pumping sum is negative
        total_periurban_repda = REPDA_array[REPDA_array[:,5]==0,4].sum() * -1 # Pumping sum is negative
        
        # Include only municipalities with urban land cover
        mun = np.unique(self.mun)[1:].copy()
        
        # Loop through model phases
        for i, l in enumerate(LU_PAR):
            # Calculate unaccounted for water supply by subtraction to determine pumping in REPDA dataset
            total_mthly_pumping = self.twateruse[i] - new_other[i] # Monthly pumping is equal to the total water use minus other supply (m3/d)
            
            LEAK_array = np.zeros((LU[l]['LIST']['URBAN'].shape[0] * (PHASE_PER[i + 1] - PHASE_PER[i]),5))
            j = 0
            
            # Generate monthly pumping datasets for REPDA data in single pumping value format
            # Loop through monthly periods in each model phase
            for p in range(S_per[i]+1,E_per[i]+1):
                unknown_pumping = total_mthly_pumping - (-1 * np.sum(list(zip(*WEL_INFO[p-1]))[3])) # Unknown pumping is the total monthly pumping for each model period minus the known pumping from SACM and CAEM (which are negative)
                                
                # Urban wells
                WEL_DICT, WEL_INFO = self.addNewWells(New_WEL=REPDA_array[REPDA_array[:,5]==1,:5], LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=self.params['Q'][i], start=p, end=(p + 1), pumpwell=True)                
                
                p_multiplier = (unknown_pumping - total_urban_repda*self.params['Q'][i]*self.ratiogn['PHASE'][i])/total_periurban_repda # Determine the monthly multiplier by dividing the estimated unknown pumping by the total pumping in the REPDA dataset
                
                # Peri-urban wells
                WEL_DICT, WEL_INFO = self.addNewWells(New_WEL=REPDA_array[REPDA_array[:,5]==0,:5], LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=p_multiplier, start=p, end=(p + 1), pumpwell=True)

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
                
                for n,m in enumerate(mun):
                    # Create an array for all urban model cells in this municipality
                    tempLeak = LU[l]['LIST']['URBAN'][(LU[l]['LIST']['URBAN'][:,4]==m),:4]
                    if len(tempLeak) > 0:
                        u_cells = tempLeak[:, 2].sum() # Number of urban model cells in municipality m
                        # Use total pumping for each stress period to determine leak quantities
                        LperCell = float(leaks[np.where(leaks==m)[0],i+1]) * (1 - fixleak) * self.params['IN'][0] / u_cells # Reference the leaks in the municipality (m3/d) multiply by (1 - fixleak), multiply by infiltration rate, divide by the number of urban cells
                        tempLeak[:,2] *= LperCell
                        
                        # apply 90% returns to sewer under clay layer (Geologic formation 1)
                        tempLeak[tempLeak[:,3]==1, 2] *= 0.1
                        
                        # Get rows of all cells of urban land use type from list
                        LEAK_array[j:(j + tempLeak.shape[0]), 0] = tempLeak[:, 0]
                        
                        # Get columns of all cells of urban land use type from list
                        LEAK_array[j:(j + tempLeak.shape[0]), 1] = tempLeak[:, 1]
                        
                        # Set the period to the current stress period for all urban cells
                        LEAK_array[j:(j + tempLeak.shape[0]), 2] = p
                        
                        # Set the end of the period to the next stress period
                        LEAK_array[j:(j + tempLeak.shape[0]), 3] = p + 1
                        
                        # Set the multiplier to the percentage of urban land use type stored in list
                        LEAK_array[j:(j + tempLeak.shape[0]), 4] = tempLeak[:,2]
                        
                        # Set the new index to the previous index plus the number of cells added
                        j += tempLeak.shape[0]
          
            WEL_DICT, WEL_INFO = self.addNewWells(LEAK_array, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=self.params['LK'][i], coordType='rc')
        
        total_mthly_leak = np.zeros(PHASE_PER[3])
        for i in range(total_mthly_leak.shape[0]):
            total_mthly_leak[i] = np.sum(list(zip(*WEL_INFO[i]))[3])
        total_mthly_leak = (total_mthly_leak - total_mthly_pumping)
        ### Cost attributed to fixing leaks
        average_leak = total_mthly_leak.sum() / (360 * sec2day) # convert to m3/s
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
        WWTPs = np.loadtxt(r'data_processed\alternatives\WWTP.csv', delimiter=',', skiprows=1, usecols=[9,8,5,6,11])
        
        if num_WWTP > 0:
            WWTPs[:, 3] = (WWTPs[:,2] - WWTPs[:,3]) # Col 2 is installed treatment capacity and Col 3 is actual treatment quantity
            WWTPs = np.insert(WWTPs, 2, PHASE_PER[0], axis=1) # Insert starting period
            WWTPs[:, 3] = np.ones(WWTPs.shape[0]) * PHASE_PER[phases] # Replace Col 3 with ending period
            WWTPs[WWTPs[:, 4] < 0.01, 4] = 0.01 # For any WWTPs with a difference of
            # less than 0.01 m3/s in capacity, assign an injection of 0.01 m3/s
            
            WWTPs = WWTPs[np.random.choice(WWTPs.shape[0], size=num_WWTP, replace=False), :]
            WEL_DICT, WEL_INFO = self.addNewWells(WWTPs, LYR=1, WEL_Dict=WEL_DICT, INFO_Dict=WEL_INFO, WEL_mult=sec2day) # Mult used to convert to m3/d
            
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
        Basin_Array = np.loadtxt(r'data_processed\alternatives\RCH_BASIN.csv', delimiter=',', skiprows=1)
        Basins = np.zeros((num_RCHBASIN, 2))
        for b in range(0, num_RCHBASIN):
            randBasin = np.random.randint(0, Basin_Array.shape[0])
            c = int(np.floor((Basin_Array[randBasin,0] - self.xll) / self.cellsize))
            r = int(np.floor((self.yur - Basin_Array[randBasin, 1]) / self.cellsize))
            Basins[b, :] = [r, c]
            
            for per in range(0, 360):
                # Add injection equal to treatment capacity for each parameter period
                WEL_DICT[per].append([1, r, c, 1 * sec2day]) # 1 m3/s recharge basins (35 cfs)
                WEL_INFO[per].append([1, r, c, 1 * sec2day, self.mun[r, c], 1])
        
        ## Cost attributed to building recharge basins
        cost += num_RCHBASIN * 20
        
        wel = flopy.modflow.ModflowWel(mf, stress_period_data=WEL_DICT, options=['NOPRINT'])
            
        print('WEL_Dict generated in', str(time.time() - newtime), 'seconds')
        
        ## Create pickle file of Well Data to be available for post processing of well energy use objective
        winfofile = r'model_files\optimization_data\objectives\WEL_INFO_' + self.name + '.pickle'
        with open(winfofile, 'wb') as handle:
            pickle.dump(WEL_INFO, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print('WEL_Dict saved in',str(time.time()-newtime),'seconds')
        
#        ## Drains
#        drain_data = np.loadtxt(r'data_processed\drains\DRN_order3.csv', delimiter=',', skiprows=1)
#        DRN_DICT = {}
#        if drains:
#            drain_list = []
#            for drn in range(0,drain_data.shape[0]):
#                r = int(drain_data[drn,0])
#                c = int(drain_data[drn,1])
#                drain_list.append([0,r,c,self.dem[r,c],drain_data[drn,2]*self.params['DRN'][0]])
#            for per in range(PHASE_PER[phases]):
#                DRN_DICT[per] = drain_list
#            self.drains = DRN_DICT
#            
#            drn = flopy.modflow.ModflowDrn(mf, stress_period_data=DRN_DICT, options=['NOPRINT'])
#                    
        # Generate output control and solver MODFLOW packages 
        oc, pcg = self.outputControl(mf)
        
        mf.add_existing_package(r'model_files\modflow\OBS.ob_hob',ptype='HOB', copy_to_model_ws=False)
        mf.add_output(r'model_files\modflow\VM_Test.hob.out',unit=1002)
        hob = flopy.modflow.ModflowHob.load(r'model_files\modflow\OBS.ob_hob', mf)
        winfofile = r'model_files\modflow\OBS.pickle'
        with open(winfofile, 'wb') as handle:
            pickle.dump(hob.obs_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
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
        success, buff = mf.run_model()
        
        print('MODFLOW model completed run in', str(time.time() - newtime), 'seconds')
        
        print('Wrapper completed run in', str(time.time() - timestart), 'seconds')
        
        self.wwtps = WWTPs
        self.basins = Basins
        self.mthlyleak = total_mthly_leak
        self.cost = cost
        self.wells = WEL_INFO
        self.landuse = LU

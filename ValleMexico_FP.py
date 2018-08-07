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

xll = 455000
yll = 2107000
xur = 539000
yur = 2175000
cellsize = 500

# Load datasets
ACTIVE = gen.openASC('data_output\ACTIVE_VM.asc')
GEO = gen.openASC('data_output\GEO_VM.asc')
DEM = gen.openASC('data_output\DEM_VM.asc')

# Assign name and create modflow model object
modelname = 'VM_Test'
#mf = flopy.modflow.Modflow.load(r'ValleyMexico_zones.nam',exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')
mf = flopy.modflow.Modflow(modelname, exe_name=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe')

# Model domain and grid definition
ztop = DEM # Model top elevation (Digital Elevation Model)
L1botm = ztop - 50 # Layer 1 bottom elevation
L2botm = L1botm - 350 # Layer 2 bottom elevation
L3botm = L2botm - 1500 # Layer 3 bottom elevation
botm = [L1botm,L2botm,L3botm] # Model bottom elevation
nlay = 3 # Number of layers
Lx = xur - xll # Model width (X length)
Ly = yur - yll # Model height (Y length)
ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns
delr = cellsize # Row height
delc = cellsize # Column height

#%% Process parameter data

#PhaseParams = np.zeros((7,))
# COL1 : Phase year
LU = [1985,1990,1995,2000,2005,2010,2015]
# COL2 : Phase period starting stress period
PHASE_PER = [1, 48, 204, 252, 324, 360]
#***** COL3 : Phase well pumping multiplier
WEL_PAR = [8.871e-1,8.626e-1,9.509e-1,8.48e-1,8.538e-1]
# COL4 : Phase distribution system leak multiplier
LEAK_Q = 60.3*86400 # m3/d
LEAK_PAR = np.array([4.28E-02, 4.85E-02, 4.86E-02, 6.00E-02, 5.75E-02])*LEAK_Q
# COL5 : Phase LID increase multiplier
LID_PAR = [1, 1, 1, 1, 1, 1, 1]
# COL : 
#ZoneParams = np.zeros((6,))
# COL1 : Zone hydraulic conductivity vertical anisotropy
VK_PAR = [1.0e2,1.0e2,1.0e-1,1.0e-1,1.0e1]
# COL2 : Zone specific storage
SS_PAR = [3.28E-03,6.56E-05,3.28E-05,3.28E-06,3.28E-07]
# COL3 : Zone hydraulic conductivity
HK_PAR = [4.32E-03,2.592E+01,4.32E+00,8.64E-02,8.6400E-05] # m/d
# COL4 : Zone specific yield
#SY_PAR = [0.06,0.15,0.15,0.01,0.01]

RCH_PAR = [3.65E-03, 3.94E-01, 1.29E-01] # Recharge multiplier for urban, natural, and irrigated cover

#%% Assign Land Use Types
URBAN = {}
NATURAL = {}
WATER = {}
for n in range(0,7):
    filename = r'data_output\LU\LU_' + str(LU[n]) + '.asc'
    for band, LUtype in enumerate([URBAN,NATURAL,WATER]):
        LU_Array = gen.openASC(filename,band)
        LUtype[LU[n]][0] = LU_Array
        
        i = 0
        node = np.zeros((nrow*ncol,7),dtype=np.int)
        node[:,0] = np.arange(0,nrow*ncol)
        
        for r in range(nrow):
            for c in range(ncol):
                node[i,1] = ACTIVE[r,c] # Determine if ACTIVE
                node[i,2] = r
                node[i,3] = c
                node[i,4] = LU_Array[r,c] # Assign LU value
                i+=1
        
        node = node[node[:,1]==1,:]
        LUtype[LU[n]][1] = node[node[:,4]==band,:]


##%% Time discretization
#nper = 360 # Number of stress periods
#nstp  = []
#for y in range(1984,2014):
#    for m in range(1,13):
#        nstp.append(calendar.monthrange(y,m)[1])
#nstp[0] = 1
#nstp = np.array(nstp)
#steady = np.zeros((360),dtype=bool)
#steady[359] = True
#
##%% Model Boundaries & initial conditions
#
## Create the discretization object
#dis = flopy.modflow.ModflowDis(mf, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, delr=delr, delc=delc,
#                               top=ztop, botm=botm, perlen = nstp, nstp=nstp, steady=steady)
#
## Active areas
#ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
#for i in range(0,3): ibound[i,:,:] = ibound[i,:,:]*ACTIVE
#ibound[0,:,:] *= (np.array(GEO==1)+np.array(GEO==2))*1
#
## Variables for the BAS package
#strt = np.zeros((nlay, nrow, ncol), dtype=np.float32)
#IH = bf.HeadFile('data\Input\IH-from-SS.hds')
#strt = IH.get_data(totim=1.0)
#
#bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt, ifrefm=True, ichflg=True, stoper=1)
#    
##%% Layer properties
## Add LPF package to the MODFLOW model
## Hydraulic conductivity
#HK_LYR1 = (GEO==1)*HK_PAR[0] + (GEO==2)*HK_PAR[1]
#HK_LYR2 = (GEO==1)*HK_PAR[1] + (GEO==2)*HK_PAR[1] + (GEO==3)*HK_PAR[2] + (GEO==4)*HK_PAR[3]
#HK_LYR3 = HK_PAR[4]
#HK = np.array([HK_LYR1,HK_LYR2,HK_LYR3])
#
## Vertical anisotropy (H:V) of hydraulic conductivity
#VK_LYR1 = (GEO==1)*VK_PAR[0] + (GEO==2)*VK_PAR[1]
#VK_LYR2 = (GEO==1)*VK_PAR[0] + (GEO==2)*VK_PAR[1] + (GEO==3)*VK_PAR[2] + (GEO==4)*VK_PAR[3]
#VK_LYR3 = VK_PAR[4]
#VKA = np.array([VK_LYR1,VK_LYR2,VK_LYR3])
#
## Specific storage
#SS_LYR1 = (GEO==1)*SS_PAR[0] + (GEO==2)*SS_PAR[1]
#SS_LYR2 = (GEO==1)*SS_PAR[0] + (GEO==2)*SS_PAR[1] + (GEO==3)*SS_PAR[2] + (GEO==4)*SS_PAR[3]
#SS_LYR3 = SS_PAR[4]
#SS = np.array([SS_LYR1,SS_LYR2,SS_LYR3])
#
### Specific yield
##SY_LYR1 = (GEO==1)*SY_PAR[0] + (GEO==2)*SY_PAR[1]
##SY_LYR2 = (GEO==1)*SY_PAR[0] + (GEO==2)*SY_PAR[1] + (GEO==3)*SY_PAR[2] + (GEO==4)*SY_PAR[3]
##SY_LYR3 = SY_PAR[4]
##SY = np.array([SY_LYR1,SY_LYR2,SY_LYR3])
#
#lpf = flopy.modflow.mflpf.ModflowLpf(mf, ipakcb=9, hdry=-1e+20, laytyp=[0,0,0], layvka=[1,1,1], 
#                                     laywet=[0,0,0], hk=HK, vka=VKA, ss=SS)
##lpf = flopy.modflow.ModflowLpf.load(r'Tlalpan1.lpf', mf)
#
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
#
##%% Plotting
#
#fig, axes = plt.subplots(2, 3, figsize=(11, 8.5))
#axes = axes.flat
#for i, hdslayer in enumerate(hds):
#    im = axes[i].imshow(hdslayer, vmin=0, vmax=hds.max())
#    axes[i].set_title('Layer {}'.format(i+1))
#    ctr = axes[i].contour(hdslayer, colors='k', linewidths=0.5)
#    
#    # export head rasters 
#    # (GeoTiff export requires the rasterio package; for ascii grids, just change the extention to *.asc)
#    mf.sr.export_array('data/heads{}.tif'.format(i+1), hdslayer)
#    
#    # export head contours to a shapefile
#    mf.sr.export_array_contours('data/heads{}.shp'.format(i+1), hdslayer)
#    
#fig.delaxes(axes[-1])
#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
#fig.colorbar(im, cax=cbar_ax, label='Head')
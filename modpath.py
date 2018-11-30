# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 20:57:53 2018

@author: MM
"""

import numpy as np
import flopy
from gwscripts.dataprocessing import gengriddata as gen
from flopy.utils.zonbud import ZoneBudget, read_zbarray

# Model domain and grid definition
xll = 25
yll = 25
xur = 525
yur = 525
cellsize = 50
ncol = int((xur-xll)/cellsize) # Number of rows
nrow = int((yur-yll)/cellsize) # Number of columns
nlay = 2 # Number of layers

m = flopy.modflow.Modflow.load('OBJ_Test.nam')
m.get_package_list()

# Porosity
GEO = gen.openASC('data_output\objTest\GEO_VM.asc')
ACTIVE1 = gen.openASC('data_output\objTest\ACTIVE_VM.asc')
N_PAR = [0.5,0.2,0.15,0.15,0.05]
p = np.zeros((nlay,nrow,ncol))
p[0,:,:] = (ACTIVE1==1)*N_PAR[0]
p[1,:,:] = (GEO==2)*N_PAR[1] + (GEO==3)*N_PAR[2] + (GEO==4)*N_PAR[3] + (GEO==5)*N_PAR[4]

zonbud_zones = np.zeros((nlay,nrow,ncol), dtype='int16')
zonbud_zones[:,0:5,0:5] = 1
zonbud_zones[:,0:5,5:10] = 2
zonbud_zones[:,5:10,0:5] = 3
zonbud_zones[:,5:10,5:10] = 4

zb = ZoneBudget(r'OBJ_Test.cbc', zonbud_zones)
df = zb.get_dataframes()
#
#mf = flopy.modflow.Modflow.load('OBJ_Test.nam')
#mp = flopy.modpath.Modpath(modelname='modpathtest',modflowmodel=mf)
#mpbas = flopy.modpath.ModpathBas(mp,prsity=p)
#mpsim = mp.create_mpsim(simtype='endpoint', trackdir='forward', packages='WEL')
#mpsim = flopy.modpath.ModpathSim(mp)
#mp.run_model()

#MODPATH
SimulateionType = 1
TrackingDirection = 1
WeakSinkOption = 1
WeakSourceOption = 1
ReferenceTimeOption = 1
StopOption = 1
ParticleGenerationOption = 1
TimePointOption = 1
BudgetOutputOption = 1
ZoneArrayOption = 1
RetardationOption = 1
AdvectiveObservationOption = 1
options = [SimulateionType,TrackingDirection,WeakSinkOption,WeakSourceOption,
           ReferenceTimeOption,StopOption,ParticleGenerationOption,TimePointOption,
           BudgetOutputOption,ZoneArrayOption,RetardationOption,AdvectiveObservationOption]

dis_file = 'OBJ_Test.dis'
head_file = 'OBJ_Test.hds'
budget_file = 'OBJ_Test.cbc'

mp = flopy.modpath.Modpath(modelname=m.name,modflowmodel=m,dis_file=dis_file,
                       head_file=head_file,budget_file=budget_file)

mnf = 'OBJ_Test.mpnam'
mlf = 'OBJ_Test.mplst'
mef = 'OBJ_Test.mpend'

mpsim = flopy.modpath.ModpathSim(mp,mp_name_file=mnf,mp_list_file=mlf,option_flags=options,extension='mpsim')

mpb = flopy.modpath.ModpathBas(mp,prsity=p)

mp.write_input()


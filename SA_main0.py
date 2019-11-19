import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
import numpy as np
#from mpi4py import MPI

# MODFLOW File location
#exefile = '~/MODFLOW/MF2005.1_12u/make/mf2005'
exefile=r'C:\WRDAPP\MF2005.1_12\bin\mf2005.exe'

# Load alternatives
altpath = Path.cwd() / 'model_files' / 'optimization_data' / 'decisions' / 'alternatives.csv'
alternatives = {}
altkeys = ['names','wwtps','basins','leakrepair','labels']
with open(altpath) as a:
    lines = a.read().splitlines() 
    for i, line in enumerate(lines):
        templist = line.split(',')
        alternatives[altkeys[i]] = templist

nsamples = 1000
soswrlim = 5000000

params = {}

# Test of JHyd parameters
#params['names'] = ['HK_Par1','HK_Par2','HK_Par3','HK_Par4','HK_Par5','SS_Par1','SS_Par2','SS_Par3','SS_Par4','SS_Par5','VANI_Par1','VANI_Par2','VANI_Par3','VANI_Par4','VANI_Par5','Q_Par1','Q_Par2','Q_Par3','LK_Par1','LK_Par2','LK_Par3','IN_Par1','TWU_Par1','TWU_Par2','TWU_Par3','RCH_Par1','RCH_Par2','RCH_Par3','SY_Par1','SY_Par2','SY_Par3','SY_Par4','SY_Par5']
#params['values'] = [[-2.99573,4.60517,-0.07085,-1.05872,-2.99573,-3.91202,-9.79159,-13.12236,-13.12236,-13.12236,1.60944,2.30259,0.00000,2.30259,2.30259,2.65300,1.46700,0.46390,1.00000,1.00000,1.00000,0.10000,0.83490,1.17500,1.20200,-4.60517,-1.02666,-0.69315,-2.81341,0.13110,-1.98195,0.26580,-4.60517]]
#parames['transform'] = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,1]

# Test from sample set
params['names'] = ['HK_Par1', 'HK_Par2', 'HK_Par3', 'HK_Par4', 'HK_Par5', 'SS_Par1', 'SS_Par2', 'SS_Par3', 'SS_Par4', 'SS_Par5', 'SY_Par1', 'SY_Par2', 'SY_Par3', 'SY_Par4', 'SY_Par5', 'VANI_Par1', 'VANI_Par2', 'VANI_Par3', 'VANI_Par4', 'VANI_Par5', 'Q_Par1', 'Q_Par2', 'Q_Par3', 'LK_Par1', 'LK_Par2', 'LK_Par3', 'TWU_Par1', 'TWU_Par2', 'TWU_Par3', 'RCH_Par1', 'RCH_Par2', 'RCH_Par3', 'IN_Par1']
params['values'] = [[-7.89307,3.30628,-1.12245,-0.679616,3.89708,-5.20895,-7.4379,-12.6959,-15.4405,-13.6696,-5.87041,0.120031,-3.54277,0.140016,-2.41413,3.2686,4.93379,1.08031,1.69464,1.74949,2.12159,2.10619,3.62973,1.51714,0.508849,1.25122,1.33426,0.835494,1.3121,-2.54135,-2.3958,-1.7944,0.300927]]
params['transform'] = [1,1,1,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,1,1,1,0]                     

# Test of SALib generated params
#params = sa.gen_param_vals(nsamples)
#np.savetxt('params_transform.csv', params['values'], delimiter=',')

# Test of csv params file from JHyd
#params = Path.cwd() / 'model_files' / 'modflow' / 'params.pval'

safolder = '0014'
sarun = 0

objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=True, delfolder=False)
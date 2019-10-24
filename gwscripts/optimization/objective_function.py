# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 21:33:40 2019

@author: MM
"""
from ValleMexico_setup import *
import flopy.utils.binaryfile as bf
from gwscripts.optimization import measureobjectives as mo
from pathlib import Path

def objective_function(x):
    num_WWTP = x[0]
    num_Basin = x[1]
    fix_leak = x[2]
    
    optpath = Path.cwd().parent.parent
    
    optmodel = model('Opt', 455000, 2107000, 539000, 2175000, 500, 1984, 2014, ACTIVE=[optpath / 'data_processed' / 'ACTIVE_VM_LYR2.asc', optpath / 'data_processed' / 'ACTIVE_VM_LYR2.asc'], THICKNESS=[optpath / 'data_processed' / 'THICK1_VM.asc', optpath / 'data_processed' / 'THICK2_VM.asc'], GEO=[optpath / 'data_processed' / 'GEO_VM_LYR1.asc', optpath / 'data_processed' / 'GEO_VM_LYR2.asc'], DEM=optpath / 'data_processed' / 'DEM_VM.asc', IH=optpath / 'data_processed' / 'IH_1984_LT2750.asc', MUN=optpath / 'data_processed' / 'MUN_VM.asc', PAR=optpath / 'model_files' / 'modflow' / 'params.pval')
    optmodel.run_simulation_model(num_WWTP,num_Basin,fix_leak)
    optmodel.heads = bf.HeadFile(optpath / 'model_output' / 'VM_Opt.hds')
    
    energy = mo.measureEnergy(optmodel.heads,optmodel.wells,optmodel.dem)
    subs_array = mo.measureSubidence(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.th1)
    mound_array = mo.measureMound(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.landuse,[132,252])
    subs = subs_array[0]/subs_array[1]
    mound = mound_array[0]
    
    return [energy,subs,mound,optmodel.cost]
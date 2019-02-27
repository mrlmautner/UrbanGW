# -*- coding: utf-8 -*-
"""
Created on Tue Feb 26 21:33:40 2019

@author: MM
"""
from ValleMexico_setup import *
import flopy.utils.binaryfile as bf
from gwscripts.optimization import measureobjectives as mo

def objective_function(x):
    num_WWTP = x[0]
    num_Basin = x[1]
    fix_leak = x[2]
    
    optmodel = model('Opt', 455000, 2107000, 539000, 2175000, 500, 1984, 2014, 'data_output\ACTIVE_VM_LYR1.asc', 'data_output\ACTIVE_VM_LYR2.asc', 'data_output\THICK1_VM.asc', 'data_output\THICK2_VM.asc', 'data_output\GEO_VM.asc', 'data_output\DEM_VM.asc', 'data_output\IH_1984.asc','data_output\MUN_VM.asc')
    optmodel.run_scenario_model(num_WWTP,num_Basin,fix_leak)
    optmodel.heads = bf.HeadFile('model_output\VM_Opt.hds')
    
    energy = mo.measureEnergy(optmodel.heads,optmodel.wells,optmodel.dem)
    subs_array = mo.measureSubidence(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.th1)
    mound_array = mo.measureMound(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.landuse,[132,252])
    subs = subs_array[0]/subs_array[1]
    mound = mound_array[0]
    
    return [energy,subs,mound,optmodel.cost]
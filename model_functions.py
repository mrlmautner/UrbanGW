# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vm
import plot_results as pltvm
from gwscripts.optimization import measureobjectives as mo
import flopy.utils.binaryfile as bf
import numpy as np
import time
import pandas as pd
from pathlib import Path
import shutil
#import pickle

# Function to run SA: sarun is the index of the set of parameter values to use given a previously generated set of parameter values in SALib, soswrlim is the maximum sum-of-squared, weighted residuals allowed to proceed with an evaluation of the managed aquifer recharge alternatives
def SA_mode(alternatives, params, exefile, safolder, sarun=0, soswrlim=0, verbose=False, delfolder=True):
    
    # Load municipality raster
    mun = np.loadtxt(Path.cwd() / 'data_processed' / 'MUN_VM.asc', skiprows=6)
    munlist = np.unique(np.unique(mun)) # List of municipalities
    munlist = munlist[munlist>0]
    numMun = munlist.shape[0]
    
    num_alts = len(alternatives['names'])
    
    # Create a folder for the sarun parameter set
    sa_loc = Path.cwd() / 'model_files' / 'modflow' / str(sarun)
    sa_loc.mkdir(exist_ok=True)
    
    # Run the simulation model under historical conditions
    hist_time = time.time()
    sa_model = vm.model(name=str(Path(sa_loc) / alternatives['names'][0]), PAR=params, sarun=sarun, exe_file=exefile)
    sa_model.run_simulation_model(alt_wwtp=0, alt_basin=0, alt_leak=0, incl_obs=True, verbose=verbose)
    print('Historical alternative for model ' + '{:05d}'.format(sarun) + ' completed in: ' + str(time.time() - hist_time) + ' seconds', flush=True)
    
    # Load head observation information for historical model run
    stats = np.loadtxt(Path.cwd() / 'model_files' / 'modflow' / 'OBS_stats.csv',skiprows=1)
    df = pd.read_fwf(Path.cwd().joinpath('model_files').joinpath('modflow').joinpath(str(sarun)).joinpath(alternatives['names'][0]+'.hob.out'),widths=[22,19,22])
    heads_obs = [df.columns.values.tolist()] + df.values.tolist()
    
    soswr, maxerror = mo.calculate_SOSWR(heads_obs, stats)
            
    # Move the head observation file into the SA experiment folder (safolder)
    sa_hob_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('hob').joinpath(safolder)
    try:
        shutil.move(str(sa_loc.joinpath(alternatives['names'][0] + '.hob.out')), str(sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out')))
    except:
        sa_hob_loc.mkdir(exist_ok=True)
        shutil.move(str(sa_loc.joinpath(alternatives['names'][0] + '.hob.out')), str(sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out')))
        
    # Save the sum-of-squared, weighted residual error
    error = np.array([soswr, maxerror])
    sa_err_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('err').joinpath(safolder)
    try:
        np.savetxt(sa_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
    except:
        sa_err_loc.mkdir(exist_ok=True)
        np.savetxt(sa_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
        
    # Determine if historical model performance is adequate     
    energy, wq, mound = (np.ones(num_alts)*np.nan for i in range(3))
    energy_MUN, wq_MUN, mound_MUN = ([[] for _ in range(num_alts)] for i in range(3))
    #heads = pltvm.get_heads([alternatives['names'][0]], sarun)
    #energy[0], wq[0], mound[0] = mo.get_objectives(heads[alternatives['names'][0]], sa_model.wells, sa_model.landuse, sa_model.dem, sa_model.actv, sa_model.botm)

    try:
        heads = pltvm.get_heads([alternatives['names'][0]], sarun)
        energy[0], energy_MUN[0], wq[0], wq_MUN[0], mound[0], mound_MUN[0] = mo.get_objectives(heads[alternatives['names'][0]], sa_model.wells, sa_model.landuse, sa_model.dem, sa_model.actv, sa_model.botm, sa_model.mun)
    except:
        energy[0], energy_MUN[0], wq[0], wq_MUN[0], mound[0], mound_MUN[0] = np.nan, np.ones(numMun)*np.nan, np.nan, np.ones(numMun)*np.nan, np.nan, np.ones(numMun)*np.nan

    # Execute the MODFLOW model for each alternative and collect results
    for i, name in enumerate(alternatives['names'][1:]):
        if verbose: print(name, 'Alternative', flush=True)
        alt_time = time.time()
        sa_model = vm.model(name=str(Path(sa_loc) / name), PAR=params, sarun=sarun, exe_file=exefile)
        sa_model.run_simulation_model(int(alternatives['wwtps'][i+1]), int(alternatives['basins'][i+1]), int(alternatives['leakrepair'][i+1]), verbose=verbose)
        if verbose: print(name, 'Simulation completed in', str(time.time() - alt_time), 'seconds', flush=True)
        print('Alternative ' + name + ' for model ' + '{:05d}'.format(sarun) + ' completed in: ' + str(time.time() - alt_time) + ' seconds', flush=True)

        # Calculate objective performance
        if verbose: print('Calculating alternative performance under objectives', flush=True)
        
        #heads = pltvm.get_heads([name], sarun)
        #energy[i+1], wq[i+1], mound[i+1] = mo.get_objectives(heads[name], sa_model.wells, sa_model.landuse, sa_model.dem, sa_model.actv, sa_model.botm)
    
        try:
            heads = pltvm.get_heads([name], sarun)
            energy[i+1], energy_MUN[i+1], wq[i+1], wq_MUN[i+1], mound[i+1], mound_MUN[i+1] = mo.get_objectives(heads[name], sa_model.wells, sa_model.landuse, sa_model.dem, sa_model.actv, sa_model.botm, sa_model.mun)
        except:
            energy[i+1], energy_MUN[i+1], wq[i+1], wq_MUN[i+1], mound[i+1], mound_MUN[i+1] = np.nan, np.ones(numMun)*np.nan, np.nan, np.ones(numMun)*np.nan, np.nan, np.ones(numMun)*np.nan
        
    mound = mound*100
    wq = wq*100
    mound_MUN = [i * 100 for i in mound_MUN]
    wq_MUN = [i * 100 for i in wq_MUN]
        
    objectives = [energy, wq, mound]
    sa_obj_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('obj').joinpath(safolder)
    try:
        np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
    except:
        sa_obj_loc.mkdir(exist_ok=True)
        np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')

    objectives_MUN = np.array([energy_MUN, wq_MUN, mound_MUN])
    sa_mun_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('mun').joinpath(safolder).joinpath(str(sarun))
    for m, mun in enumerate(munlist):
        try:
            np.savetxt(sa_mun_loc.joinpath('{:05.0f}'.format(mun) + '.csv'), objectives_MUN[:,:,m], delimiter=',')
        except:
            sa_mun_loc.mkdir(parents=True,exist_ok=True)
            np.savetxt(sa_mun_loc.joinpath('{:05.0f}'.format(mun) + '.csv'), objectives_MUN[:,:,m], delimiter=',')

#    # Save last year of heads
#    sa_heads_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('heads').joinpath(safolder)
#    try:
#        f = open(sa_heads_loc.joinpath('{:05d}'.format(sarun) + '.pkl'),"wb")
#        pickle.dump(h,f)
#        f.close()
#    except:
#        sa_heads_loc.mkdir(exist_ok=True)
#        f = open(sa_heads_loc.joinpath('{:05d}'.format(sarun) + '.pkl'),"wb")
#        pickle.dump(h,f)
#        f.close()

    # Delete the model directory and files for this SA parameter set
    if delfolder:
        shutil.rmtree(str(sa_loc), ignore_errors=True)

    return [error, objectives]

#def checkHOB():
#    Path.cwd().joinpath('Input').joinpath('params').joinpath('params_' + safolder + '.csv')
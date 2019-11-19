# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vm
import plot_results as pltvm
#from gwscripts.optimization import opttools as opt
from gwscripts.optimization import measureobjectives as mo
import flopy.utils.binaryfile as bf
import numpy as np
import time
#import pickle
import pandas as pd
from pathlib import Path
#from platypus import Problem, Integer, Real, NSGAII
#import os
import shutil

# Function to run SA: sarun is the index of the set of parameter values to use given a previously generated set of parameter values in SALib, soswrlim is the maximum sum-of-squared, weighted residuals allowed to proceed with an evaluation of the managed aquifer recharge alternatives
def SA_mode(alternatives, params, exefile, safolder, sarun=0, soswrlim=0, verbose=False, delfolder=True):
    num_alts = len(alternatives['names'])
    sa_model = [0]*num_alts
    
    # Create a folder for the sarun parameter set
    sa_loc = Path.cwd() / 'model_files' / 'modflow' / str(sarun)
    sa_loc.mkdir(exist_ok=True)
    
    # Run the simulation model under historical conditions
    hist_time = time.time()
    sa_model[0] = vm.model(name=str(Path(sa_loc) / alternatives['names'][0]), PAR=params, sarun=sarun, exe_file=exefile)
    sa_model[0].run_simulation_model(alt_wwtp=0, alt_basin=0, alt_leak=0, incl_obs=True, verbose=verbose)
    print('Historical alternative for model ' + '{:05d}'.format(sarun) + ' completed in: ' + str(time.time() - hist_time) + ' seconds')
    
    # Load head observation information for historical model run
    stats = np.loadtxt(Path.cwd() / 'model_files' / 'modflow' / 'OBS_stats.csv')
    df = pd.read_fwf(Path.cwd().joinpath('model_files').joinpath('modflow').joinpath(sa_loc).joinpath(alternatives['names'][0]+'.hob.out'),widths=[22,19,22])
    heads_obs = [df.columns.values.tolist()] + df.values.tolist()
    
    soswr, maxerror = mo.calculate_SOSWR(heads_obs, stats)
    
    # Save the sum-of-squared, weighted residual error
    error = np.array([soswr, maxerror])
    sa_err_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('err').joinpath(safolder)
    try:
        np.savetxt(sa_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
    except:
        sa_err_loc.mkdir(exist_ok=True)
        np.savetxt(sa_err_loc.joinpath('{:05d}'.format(sarun) + '.csv'), error, delimiter=',')
    
    # Determine if historical model performance is adequate
    objectives = np.ones((3,4))*np.nan
    if soswr<=soswrlim:
        # Execute the MODFLOW model for each alternative and collect results
        for i, name in enumerate(alternatives['names'][1:]):
            if verbose: print(name, 'Alternative')
            alt_time = time.time()
            sa_model[i+1] = vm.model(name=str(Path(sa_loc) / name), PAR=params, sarun=sarun, exe_file=exefile)
            sa_model[i+1].run_simulation_model(int(alternatives['wwtps'][i+1]), int(alternatives['basins'][i+1]), int(alternatives['leakrepair'][i+1]), verbose=verbose)
            if verbose: print(name, 'Simulation completed in', str(time.time() - alt_time), 'seconds')
            print('Alternative ' + name + ' for model ' + '{:05d}'.format(sarun) + ' completed in: ' + str(time.time() - alt_time) + ' seconds')
    
        # Retrieve head changes over model period  
        heads = pltvm.get_heads(alternatives['names'],sarun)
        
        # Calculate objective performance
        if verbose: print('Calculating alternative performance under objectives')
        energy, subs, mound = (np.zeros(num_alts) for i in range(3))
    
        for i, name in enumerate(alternatives['names']):
            energy[i], subs[i], mound[i] = mo.get_objectives(heads[name], sa_model[i].wells, sa_model[i].landuse, sa_model[i].dem, sa_model[i].geo[0], sa_model[i].thck[0], sa_model[i].botm)
        
        mound = mound*100
        
        objectives = [energy, subs, mound]
        sa_obj_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('obj').joinpath(safolder)
        try:
            np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
        except:
            sa_obj_loc.mkdir(exist_ok=True)
            np.savetxt(sa_obj_loc.joinpath('{:05d}'.format(sarun) + '.csv'), objectives, delimiter=',')
        
    # Move the head observation file into the SA experiment folder (safolder)
    sa_hob_loc = Path.cwd().joinpath('model_files').joinpath('output').joinpath('sa').joinpath('hob').joinpath(safolder)
    try:
        shutil.move(sa_loc.joinpath(alternatives['names'][0] + '.hob.out'), sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out'))
    except:
        sa_hob_loc.mkdir(exist_ok=True)
        shutil.move(sa_loc.joinpath(alternatives['names'][0] + '.hob.out'), sa_hob_loc.joinpath('{:05d}'.format(sarun) + '.hob_out'))
        
    # Delete the model directory and files for this SA parameter set
    if delfolder:
        shutil.rmtree(sa_loc, ignore_errors=True)

    return [error, objectives]

def plt_hydrograph(name, hydrographloc, exefile, verbose=True):
    hydModel = vm.model(name=name, exe_file=exefile)
    hydModel.run_simulation_model(0,0,0,verbose=verbose)
    
    pltvm.plt_wellhydrographs(name, hydrographloc, df=0, obsformation=0)
    
    return hydModel

def plt_alt(alternatives, exefile, run_alternatives=True, verbose=True):
    '''
    Alternatives mode allows the user to compare predetermined alternatives defined above by the quantity of each recharge intervention
    '''
    num_alts = len(alternatives['names'])
    alt_model = [0]*num_alts
    
    if run_alternatives:
        # Execute the MODFLOW model for each alternative and collect results
        for i, name in enumerate(alternatives['names']):
            if verbose: print(name, 'Alternative')
            alt_time = time.time()
            alt_model[i] = vm.model(name=name, exe_file=exefile)
            alt_model[i].run_simulation_model(int(alternatives['wwtps'][i]), int(alternatives['basins'][i]), int(alternatives['leakrepair'][i]), verbose=verbose)
            if verbose: print(name, 'Simulation completed in', str(time.time() - alt_time), 'seconds')
#    else:
#        # If the results have already been generated, open results from saved files
#        for i, name in enumerate(alternatives['names']):
#            if verbose: print('Opening', name, 'Alternative')
#            
#            alt_model[i] = vm.model(name=name, exe_file=exefile)
#            
#            with open(Path.cwd().joinpath('model_files').joinpath('optimization_data').joinpath('objectives').joinpath('WEL_INFO_'+name+'.pickle'), 'rb') as handle:
#                alt_model[i].wells = pickle.load(handle)
#            with open(Path.cwd().joinpath('model_files').joinpath('optimization_data').joinpath('objectives').joinpath('LU_'+name+'.pickle'), 'rb') as handle:
#                alt_model[i].landuse = pickle.load(handle)

    # Retrieve and plot head changes over model period  
    heads = pltvm.get_heads(alternatives['names'])
    
    pltvm.plt_head_change(alternatives['names'], alternatives['labels'], heads, alt_model[0].geo[1], alt_model[0].actv[1])
    
    # Calculate and plot objective performance
    if verbose: print('Calculating alternative performance under objectives')
    energy, subs, mound = (np.zeros(num_alts) for i in range(3))
    
    for i, name in enumerate(alternatives['names']):
        energy[i], subs[i], mound[i] = mo.get_objectives(heads[name], alt_model[i].wells, alt_model[i].landuse, alt_model[i].dem, alt_model[i].geo[0], alt_model[i].thck[0])
    
    mound = mound*100# mound/min(mound)
    
    pltvm.plt_alt_objectives(alternatives['names'], num_alts, [energy, subs, mound])
    
    return alt_model, energy, subs, mound
    
def objective_function(x):
    
    num_WWTP = x[0]
    num_Basin = x[1]
    fix_leak = x[2]
    
    optmodel = vm.model(name='Opt')
    optmodel.run_simulation_model(num_WWTP,num_Basin,fix_leak)
    optmodel.heads = bf.HeadFile(Path.cwd() / 'model_output' / 'VM_Opt.hds')
    
    energy = mo.measureEnergy(optmodel.heads,optmodel.wells,optmodel.dem)
    subs_array = mo.measureSubidence(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.th1)
    mound_array = mo.measureMound(optmodel.heads,optmodel.dem,optmodel.actv1,optmodel.landuse,[132,252])
    subs = subs_array[0]/subs_array[1]
    mound = mound_array[0]
    
    return [energy,subs,mound,optmodel.cost]

#def plt_opt(opt_run, max_nfes, numdec, numobj, exefile, run_optimization=True):
#    ''' 
#    The optimize option runs an optimization problem with the model using the recharge decisions and objectives defined above
#    If the run_optimization is not activated, previous results are loaded from the file indicated above
#    '''
#    if run_optimization:
#        
#        problem = Problem(numdec, numobj)
#        problem.directions[:] = Problem.MINIMIZE
#        wwtp_int = Integer(0,74) # Number of wastewater treatment plants
#        basin_int = Integer(0,10) # Number of recharge basins
#        leak_int = Integer(0,100) # 0 to 100% leak repair
#        problem.types[:] = [wwtp_int, basin_int, leak_int]
#        problem.function = objective_function()
#        algorithm = NSGAII(problem)
#        algorithm.run(max_nfes)
#    
#        results_list = []
#        variable_list = []
#        int_list = [wwtp_int, basin_int, leak_int]
#        
#        for r, results in enumerate(algorithm.result):
#            results_list.append(results.objectives[:])
#            var_list = []
#            for v,var in enumerate(results.variables):
#                var_list.append(int_list[v].decode(var))
#            variable_list.append(var_list)
#        
#        pickle.dump(results_list, open(Path.cwd().joinpath('model_files').joinpath('output').joinpath('opt').joinpath('objectives_' + opt_run + '.pkl'), "wb" ))
#        pickle.dump(variable_list, open(Path.cwd().joinpath('model_files').joinpath('output').joinpath('opt').joinpath('dvariables_' + opt_run + '.pkl'), "wb" ))
#    else:
#        with open(Path.cwd().joinpath('model_files').joinpath('output').joinpath('opt').joinpath('objectives_' + opt_run + '.pkl'), 'rb') as handle:
#            results_list = pickle.load(handle)
#        with open(Path.cwd().joinpath('model_files').joinpath('output').joinpath('opt').joinpath('dvariables_' + opt_run + '.pkl'), 'rb') as handle:
#            variable_list = pickle.load(handle)
#    
#    nondom_VM = opt.nondom_sort(results_list)
#    npresults = np.array(results_list)
#    nondom_results = npresults[nondom_VM]
#    
#    pltvm.parallel_axis(nondom_results, obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'], opt_run = opt_run)
#    
#    return results_list, variable_list, nondom_VM, npresults, nondom_results

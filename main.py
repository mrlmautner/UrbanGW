# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

from ValleMexico_setup import *
import plot_results as pltvm
from gwscripts.optimization import opttools as opt
from platypus import Problem, Integer, Real, NSGAII
import numpy as np
import pickle

# Scenario Mode
run_scenarios = False
plt_scen = False
scenario_names = ['Historical','WWTP','Leak','Basin']
leak_repair = [0,0,0.2,0]
num_wwplants = [0,74,0,0]
num_infbasins = [0,0,0,5]

# Optimization parameters
run_optimization = False
plot_opt = False
recharge_decisions = 3
recharge_objectives = 4
max_nfes = 200
opt_run = str(max_nfes)+'nfe'

'''
The model run will output the following variables:
    
w is the list of randomly selected wastewater treatment plants where reuse has
been implemented
b is the list of the row and column where each infiltration basin has been
implemented
l is an array with the total quantity of leaks in the system per month in m3
cost is the total number of interventions time their weights defined in the
model set-up
WEL_INFO is the dictionary of well objects input into the MODFLOW model which
includes positive flow from wastewater treatment plants, leaks, and recharge
basins and negative flow from pumping wells
LU is the dictionary that contains a raster and list of the percentage of land
use type per model cell for each model phase    
'''
testModel = model(455000, 2107000, 539000, 2175000, 500, 1984, 2014, 'data_output\ACTIVE_VM_LYR1.asc', 'data_output\ACTIVE_VM_LYR2.asc', 'data_output\THICK1_VM.asc', 'data_output\THICK2_VM.asc', 'data_output\GEO_VM.asc', 'data_output\DEM_VM.asc', 'data_output\IH_1984.asc','data_output\MUN_VM.asc')
wwtp, basin, leak, cost, well_info, land_use = testModel.run_scenario_model('Test',0,0,0)

if plt_scen:
    '''
    Scenario mode allows the user to compare predetermined scenarios defined
    above by the quantity of each recharge intervention
    '''
    num_scenarios = len(scenario_names)
    wwtp, basin, leak, cost, well_info, land_use = ([0]*num_scenarios for i in range(6))
    
    if run_scenarios:
        # Execute the MODFLOW model for each scenario and collect results
        for i, s_name in enumerate(scenario_names):
            vmmodel
            wwtp[i], basin[i], leak[i], cost[i], well_info[i], land_use[i] = vmmodel.run_scenario_model(s_name,
                                                             num_wwplants[i],num_infbasins[i],leak_repair[i])

    else:
        # If the results have already been generated, open results from saved files
        for i, s_name in enumerate(scenario_names):
            
            with open('model_output\objective_data\WEL_INFO_'+s_name+'.pickle', 'rb') as handle:
                well_info[i] = pickle.load(handle)
            with open('model_output\objective_data\LU_'+s_name+'.pickle', 'rb') as handle:
                land_use[i] = pickle.load(handle)
            
    s_heads = pltvm.get_heads(scenario_names)
    ### Generalize
    heads = S_heads[s_name]
    
    energy_array[s,:] = mo.measureEnergy(heads,well_info,DEM)
    subs_array[s,:] = mo.measureSubidence(heads,DEM,ACTIVE_LYR1,TH1)
    mound_array[s,:] = mo.measureMound(heads,DEM,ACTIVE_LYR1,land_use,[132,252])

if plot_opt:
    
    ''' The optimize option runs an optimization problem with the model using the
    recharge decisions and objectives defined above
    If the run_optimization is not activated, previous results are loaded from
    the file indicated above
    '''
    if run_optimization:
        
        
        nfe = 0
        problem = Problem(recharge_decisions, recharge_objectives)
        problem.directions[:] = Problem.MINIMIZE
        wwtp_int = Integer(0,74) # Number of wastewater treatment plants
        basin_int = Integer(0,10) # Number of recharge basins
        leak_int = Integer(0,100) # 0 to 100% leak repair
        problem.types[:] = [wwtp_int, basin_int, leak_int]
        problem.function = vmmodel.objective_function()
        algorithm = NSGAII(problem)
        algorithm.run(max_nfes)
    
        results_list = []
        variable_list = []
        int_list = [wwtp_int, basin_int, leak_int]
        
        for r, results in enumerate(algorithm.result):
            results_list.append(results.objectives[:])
            var_list = []
            for v,var in enumerate(results.variables):
                var_list.append(int_list[v].decode(var))
            variable_list.append(var_list)
        
        pickle.dump(results_list, open(r'model_output\opt\objectives_' + opt_run + '.pkl', "wb" ))
        pickle.dump(variable_list, open(r'model_output\opt\dvariables_' + opt_run + '.pkl', "wb" ))
    else:
        with open(r'model_output\opt\objectives_' + opt_run + '.pkl', 'rb') as handle:
            results_list = pickle.load(handle)
        with open(r'model_output\opt\dvariables_' + opt_run + '.pkl', 'rb') as handle:
            variable_list = pickle.load(handle)
    
    nondom_VM = opt.nondom_sort(results_list)
    npresults = np.array(results_list)
    nondom_results = npresults[nondom_VM]
    
    pltvm.parallel_axis(nondom_results,filename = 'parallelaxis_' + opt_run + '.png',
                        obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])


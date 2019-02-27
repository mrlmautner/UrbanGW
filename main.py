# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

from ValleMexico_setup import *
import plot_results as pltvm
from objective_function import *
from gwscripts.optimization import opttools as opt
from gwscripts.optimization import measureobjectives as mo
from platypus import Problem, Integer, Real, NSGAII
import numpy as np
import time
import pickle

# Test Mode
test = False

# Scenario Mode
plt_scen = True
run_scenarios = True
scenario_names = ['Historical','WWTP','Leak','Basin']
mapTitles = ['Historical','Increased WW Reuse','Repair Leaks','Recharge Basins']
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

if test:
    testModel = model('Test', 455000, 2107000, 539000, 2175000, 500, 1984, 2014, 'data_output\ACTIVE_VM_LYR1.asc', 'data_output\ACTIVE_VM_LYR2.asc', 'data_output\THICK1_VM.asc', 'data_output\THICK2_VM.asc', 'data_output\GEO_VM.asc', 'data_output\DEM_VM.asc', 'data_output\IH_1984.asc','data_output\MUN_VM.asc')
    testModel.run_scenario_model(0,0,0)


if plt_scen:
    '''
    Scenario mode allows the user to compare predetermined scenarios defined
    above by the quantity of each recharge intervention
    '''
    num_scenarios = len(scenario_names)
    vmmodel = [0]*num_scenarios
    
    if run_scenarios:
        # Execute the MODFLOW model for each scenario and collect results
        for i, s_name in enumerate(scenario_names):
            print(s_name + 'Scenario')
            scen_time = time.time()
            vmmodel[i] = model(s_name, 455000, 2107000, 539000, 2175000, 500, 1984, 2014, 'data_output\ACTIVE_VM_LYR1.asc', 'data_output\ACTIVE_VM_LYR2.asc', 'data_output\THICK1_VM.asc', 'data_output\THICK2_VM.asc', 'data_output\GEO_VM.asc', 'data_output\DEM_VM.asc', 'data_output\IH_1984.asc','data_output\MUN_VM.asc')
            vmmodel[i].run_scenario_model(num_wwplants[i], num_infbasins[i], leak_repair[i])
            print('Scenario' + s_name + 'completed in' + str(time.time() - scen_time) + 'seconds')
    else:
        # If the results have already been generated, open results from saved files
        for i, s_name in enumerate(scenario_names):
            vmmodel[i] = model(s_name, 455000, 2107000, 539000, 2175000, 500, 1984, 2014, 'data_output\ACTIVE_VM_LYR1.asc', 'data_output\ACTIVE_VM_LYR2.asc', 'data_output\THICK1_VM.asc', 'data_output\THICK2_VM.asc', 'data_output\GEO_VM.asc', 'data_output\DEM_VM.asc', 'data_output\IH_1984.asc','data_output\MUN_VM.asc')
            with open('model_output\objective_data\WEL_INFO_'+s_name+'.pickle', 'rb') as handle:
                vmmodel[i].wells = pickle.load(handle)
            with open('model_output\objective_data\LU_'+s_name+'.pickle', 'rb') as handle:
                vmmodel[i].landuse = pickle.load(handle)

    # Retrieve and plot head changes over model period  
    s_heads = pltvm.get_heads(scenario_names)
    
    pltvm.plt_head_change(scenario_names, mapTitles, s_heads, vmmodel[0].geo, vmmodel[0].actv2)
    
    # Calculate and plot objective performance
    energy, subs, mound = (np.zeros(num_scenarios) for i in range(3))
    
    for i, s_name in enumerate(scenario_names):
        energy[i], subs[i], mound[i] = mo.get_objectives(s_heads[i], vmmodel[i].wells, vmmodel[i].landuse, vmmodel[i].dem, vmmodel[i].actv1, vmmodel[i].th1)
    
    mound = mound/min(mound)
    
    pltvm.plt_scen_objectives(scenario_names, num_scenarios, [energy, subs, mound])
    

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
        problem.function = objective_function()
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
    
    pltvm.parallel_axis(nondom_results, filename = 'parallelaxis_' + opt_run + '.png', obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])


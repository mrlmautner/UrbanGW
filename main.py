# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vmmodel
import plot_results as pltvm
from gwscripts.optimization import opttools as opt
from platypus import Problem, Integer, Real, NSGAII
import numpy as np
import pickle

# Scenario Mode
scenario_mode = True

# Optimization parameters
optimize = False
plot_opt = True
recharge_decisions = 3
recharge_objectives = 4
max_nfes = 200

WWTPs, Basins, total_pump_leak = vmmodel.run_scenario_model('Test',0,0,100)

if scenario_mode:
    # AGU Model Runs
    numscenarios = 4
    scenarioList = ['Historical','WWTP','Leak','Basin'] 
    fixleak = [1,1,0.8,1]
    num_WWTP = [0,74,0,0]
    num_RCHBASIN = [0,0,0,5]
    w = [0]*numscenarios
    b = [0]*numscenarios
    l = [0]*numscenarios
    
    for i in range(numscenarios):
        w[i], b[i], l[i], cost, WEL_INFO, LU = vmmodel.run_scenario_model(scenarioList[i],num_WWTP[i],num_RCHBASIN[i],fixleak[i])


if optimize:
    ''' The optimize option runs an optimization problem with the model using the
    recharge decisions and objectives defined above
    '''
    
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
    
    pickle.dump(results_list, open(r'model_output\opt\objectives_200nfe.pkl', "wb" ))
    pickle.dump(variable_list, open(r'model_output\opt\dvariables_200nfe.pkl', "wb" ))
else:
    with open(r'model_output\opt\objectives_200nfe.pkl', 'rb') as handle:
        results_list = pickle.load(handle)
    with open(r'model_output\opt\dvariables_200nfe.pkl', 'rb') as handle:
        variable_list = pickle.load(handle)

if plot_opt:
    nondom_VM = opt.nondom_sort(results_list)
    npresults = np.array(results_list)
    nondom_results = npresults[nondom_VM]
    
    pltvm.parallel_axis(nondom_results,filename = 'parallelaxis'+str(max_nfes)+'.png',
                        obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])

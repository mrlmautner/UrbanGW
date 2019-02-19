# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vmmodel
import plot_results as pltvm
from gwscripts.optimization import opttools as opt
from platypus import Problem, Integer, Real, NSGAII


# Optimization parameters
optimize = False
recharge_decisions = 3
recharge_objectives = 4
max_nfes = 200

WWTPs, Basins, total_pump_leak = vmmodel.run_scenario_model('Test',0,0,100)

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

    first_variable = algorithm.result[0].variables[0]
    second_variable = algorithm.result[0].variables[1]
    third_variable = algorithm.result[0].variables[2]
    
    print(wwtp_int.decode(first_variable))
    print(basin_int.decode(second_variable))
    print(leak_int.decode(third_variable))
    
    #%%
    print(algorithm.result.variables)
    
    #%%
    results_list = []
    variable_list = []
    int_list = [wwtp_int, basin_int, leak_int]
    
    for r, results in enumerate(algorithm.result):
        results_list.append(results.objectives[:])
        var_list = []
        for v,var in enumerate(results.variables):
            var_list.append(int_list[v].decode(var))
        variable_list.append(var_list)
    
    pltvm.parallel_axis(nondom_results,filename = 'parallelaxis_200.png',
                        obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])
    
    pickle.dump(results_list, open(r'model_output\opt\objectives_200nfe.pkl', "wb" ))
    pickle.dump(variable_list, open(r'model_output\opt\dvariables_200nfe.pkl', "wb" ))

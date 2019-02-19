# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 20:27:21 2019

@author: MM
"""

import ValleMexico_setup as vmmodel
import plot_results as pltvm
from gwscripts.optimization import opttools as opt
from platypus import Problem, Integer, Real, NSGAII

WWTPs, Basins, total_pump_leak = vmmodel.run_scenario_model('Test',0,0,1)

#numscenarios = 4
#scenarioList = ['WWTP','Historical','Leak','Basin'] 
#fixleak = [1,1,0.8,1]
#num_WWTP = [74,0,0,0] # 74
#num_RCHBASIN = [0,0,0,5] # 5
#w = [0]*numscenarios
#b = [0]*numscenarios
#l = [0]*numscenarios
#
#for i in range(numscenarios):
#    w[i],b[i],l[i] = run_scenario_model(scenarioList[i],num_WWTP[i],num_RCHBASIN[i],fixleak[i])

pltvm.parallel_axis(nondom_results,filename = 'parallelaxis_200.png',
                    obj_labels = ['Energy Use\n(kWh)','Subsidence Avoidance\n(mbgs)','Urban Mounding\n(m)', 'Cost'])
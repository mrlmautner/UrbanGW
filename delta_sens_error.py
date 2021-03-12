# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 23:00:42 2021

@author: MM
"""
from pathlib import Path
import numpy as np
from SALib.analyze import delta
import gwscripts.sensitivityanalysis.satools as sa

sens_dict = {}
threshold = '0.93836'
samples = 93836
for c in ['C-00001','C-00002','C-00003','C-00004','C-00005','C-12345']:
    print(c)
    data = np.loadtxt(str(Path.cwd().joinpath('error_sets').joinpath(threshold).joinpath(c+'.csv')), delimiter=',')
    Y = data[:,33]
    
    # Generate problem set for parameters
    params = sa.gen_param_vals(samples)
    params['values'] = data[:,:33]

    sens_dict[c] = delta.analyze(params['problem'], params['values'], Y)
    
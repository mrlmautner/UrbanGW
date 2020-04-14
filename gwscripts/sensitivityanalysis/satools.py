# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:22:48 2019

@author: MM
"""

from SALib.sample import latin
from SALib.analyze import delta
import numpy as np
from pathlib import Path

def gen_param_vals(nsamples):
    prangepath = Path.cwd() / 'data_processed' / 'pranges.csv'
    # Initialize dictionary to hold model parameters
    params = {}
    paramkeys = ['names','lbound','ubound','transform']
    
    # Open text file containing parameter characteristics
    # Line 1: Parameter names
    # Line 2: Lower bound for parameter values
    # Line 3: Upper bound for parameter values
    # Line 4: Binary indicator if parameter values should be transformed (0 or 1)
    with open(prangepath) as d:
        lines = d.read().splitlines()
        for i, line in enumerate(lines):
            templist = line.split(',')
            params[paramkeys[i]] = templist
    numparams = len(params['names'])
    
    # Format lower and upper parameter bounds and apply ln transformation for applicable variables
    bounds = []
    for i, transform in enumerate(params['transform']):
        lb = float(params['lbound'][i])
        ub = float(params['ubound'][i])
        if int(transform) == 0:
            bounds.append([lb,ub])
        else:
            bounds.append([np.log(lb),np.log(ub)])   
 
    # Define the SA problem
    problem = {
        'num_vars': numparams,
        'names': params['names'],
        'bounds': bounds
    }
    
    print(problem)
    
    # Sample from the parameter ranges
    params['values'] = latin.sample(problem, nsamples)
    
    return params
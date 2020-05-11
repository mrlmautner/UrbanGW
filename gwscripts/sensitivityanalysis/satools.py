# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 11:22:48 2019

@author: MM
"""

from SALib.sample import latin
from SALib.analyze import delta
import numpy as np
from pathlib import Path
import numpy as np
import flopy
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
from scipy import stats
from pathlib import Path
import os
from SALib.sample import latin
from SALib.analyze import delta
from sklearn.cluster import KMeans

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
    
    params['problem'] = problem
    print(problem)
    
    # Sample from the parameter ranges
    params['values'] = latin.sample(problem, nsamples)
    
    return params

def kmeans_obs(obsdata, n_clusters=10, norm=False, time=True):
    '''
    obsdata is a dataframe with columns: IDPOZO, LAT, LON, Z, t
    num_cluster is the desired number of k-means clusters to group the data into
    '''
    if time:
        df_keys = ['LAT','LON','Z','t']
    else:
        df_keys = ['LAT','LON','Z']
    df = obsdata[df_keys].copy()
    if norm:
        for i in df_keys:
            df[i] = (df[i]-df[i].min())/(df[i].max()-df[i].min())
    kmeans_labels = KMeans(n_clusters=n_clusters).fit_predict(df)
    
    return kmeans_labels
    
def get_err_data(safolder, startpath, makedict=True, err_df=[], soswrlim=5000000):
    '''
    makedict option creates three dictionaries to hold error data
    otherwise, error data is returned as a dataframe
    '''
    errpath = startpath.joinpath('model_files').joinpath('output').joinpath('sa').joinpath('err')
    
    if makedict:
        # Create two arrays: soswr and max error
        soswr = {}
        maxerror = {}
        lowerror = {}
        index = 0
        
        for i in safolder:
            # Set cwd to repository folder
            errsamples = os.listdir(errpath.joinpath(i))
            numerr = len(errsamples)
            
            for x in range(numerr):
                dataFileName = errpath.joinpath(i).joinpath(errsamples[x])
                xname = int(errsamples[x][:errsamples[x].find('.csv')]) + index
                xname = '{:05d}'.format(xname)
                
                temperror = np.loadtxt(dataFileName, delimiter=',')
                soswr[xname] = temperror[0]
                
                if temperror[0] <= soswrlim:
                    lowerror[xname] = temperror[0]
                    
                maxerror[xname] = temperror[1]
                
            index += numerr # For combo sets add the number of samples thus far
            
        return soswr, maxerror, lowerror
    
    else:
        toterr = 0

        for i in safolder:
            # Set cwd to repository folder
            toterrsamples = errpath.joinpath(i)
            toterr += len(toterrsamples)
        
        err_df['soswr'] = np.ones(toterr)*np.nan
        err_df['maxerror'] = np.ones(toterr)*np.nan
        err_df['lowerror'] = np.zeros(toterr)
        index = 0
        
        for i in safolder:
            print('Processing folder:',i)
            # Set cwd to repository folder
            errsamples = os.listdir(errpath.joinpath(i))
            numerr = len(errsamples)
            
            
            for x in range(numerr):
                dataFileName = errpath.joinpath(i).joinpath(errsamples[x])
                xname = int(errsamples[x][:errsamples[x].find('.csv')]) + index
                xname = '{:05d}'.format(xname)
                
                temperror = np.loadtxt(dataFileName, delimiter=',')
                
                err_df.at[xname,'soswr'] = temperror[0]
                err_df.at[xname,'maxerror'] = temperror[1]
                
                if temperror[0] <= soswrlim:
                    err_df.at[xname,'lowerror'] = 1
                    
            index += numerr # For combo sets add the number of samples thus far
            
        return err_df

def get_obj_data(safolder, startpath, paramsarray, soswrlim=5000000, objectives=3, alternatives=['Historical', 'WWTP', 'Basin', 'Leak'], obj_names=['Energy','Quality','Flood','SOSWR','Params']):
    '''
    safolder is the folder path for the objectives to be loaded
    soswr is a dictionary with keys that correspond to each sample and contains the sum of squared weighted residual
    params array is the full dataset of parameter values for each sample
    '''
    objpath = startpath.joinpath('model_files').joinpath('output').joinpath('sa').joinpath('obj')
    
    obj = {}
    obj['arrays'] = {}
    obj['samples'] = {}
    
    # Initialize lists containing objective values
    for a in alternatives:
        obj['arrays'][a] = {}
        for o in obj_names[:objectives]:
            obj['arrays'][a][o] = []
    obj['arrays']['SOSWR'] = []
    obj['arrays']['Params'] = []
    
    index = 0
    
    for s in safolder:
        # Set cwd to repository folder
        objsamples = os.listdir(objpath / s)
        numobjsamples = len(objsamples)
        
        # Iterate over all samples that were evaluated using the planning alternatives
        for x in range(numobjsamples):
            dataFileName = objpath / s / objsamples[x]
            xname = int(objsamples[x][:objsamples[x].find('.csv')]) + index
            xname = '{:05d}'.format(xname)
            
            if soswr[xname] <= soswrlim:
                tempobj = np.loadtxt(dataFileName, delimiter=',')
                obj['samples'][xname] = tempobj
            
                # Exclude any samples that produced nan for any objectives
                if np.logical_not(np.isnan(tempobj).any()):
                    obj['arrays']['Params'].append(paramsarray[int(xname)])
                    obj['arrays']['SOSWR'].append(soswr[xname])
                    
                    for i, a in enumerate(alternatives):    
                        for j, o in enumerate(obj_names[:objectives]):
                            obj['arrays'][a][o].append(tempobj[j,i])
        
        numrunsamples = int(s[(s.find('_')+1):])
        index += numrunsamples
    
    obj['delta'] = {}
    obj['delta']['SOSWR'] = delta.analyze(problem, np.array(paramsarray), np.array(err))
    
    problem, problemtransform, params = sat.gen_param_vals(len(obj['arrays']['Params']), prangefile='pranges_full')
    
    for i, a in enumerate(altnames):
        print('Alternative:',a)
        obj['delta'][a] = {}
        for o in objnames:
            print('Objective:',o)
            obj['delta'][a][o] = delta.analyze(problem, np.array(obj['arrays']['Params']), np.array(obj['arrays'][a][o]))
            
    # Calculate Delta Sensitivity for Error
    obj['delta']['SOSWR_obj'] = delta.analyze(problem, np.array(obj['arrays']['Params']), np.array(obj['arrays']['SOSWR']))
    
    # Normalize Datasets
    obj['delta_norm'] = {}
    for i, a in enumerate(altnames):
        obj['delta_norm'][a] = {}
        for o in objnames:
            obj['delta_norm'][a][o] = obj['delta'][a][o]['delta'] / obj['delta'][a][o]['delta'].sum()
    
    obj['delta_hist_diff'] = {}
    for i, a in enumerate(altnames):
        obj['delta_hist_diff'][a] = {}
        for o in objnames:
            obj['delta_hist_diff'][a][o] = (obj['delta_norm'][a][o] - obj['delta_norm']['Historical'][o])#/obj['delta_norm']['Historical'][o]
            
    return obj

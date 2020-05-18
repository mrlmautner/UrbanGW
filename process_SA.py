# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 14:25:09 2020

@author: MM
"""

import plot_results as pltvm
import gwscripts.sensitivityanalysis.satools as sa
import gwscripts.optimization.measureobjectives as mo 
from pathlib import Path
import numpy as np
import os
import pandas as pd
from mpi4py import MPI
import itertools
from SALib.analyze import delta

folder = '20200403_100000'
startpath = Path.cwd() #Path('D:\MMautner\Cluster') #
outputpath = startpath.joinpath('model_files').joinpath('output').joinpath('sa')
hobspath = outputpath.joinpath('hob').joinpath(folder)
objpath = outputpath.joinpath('obj').joinpath(folder)

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_proc = comm.Get_size()

n_clusters = 6
n_samples = 100000
samples_per_proc = int(n_samples / num_proc)

kmeans_time = False # Include time as a variable for K-means clustering
kmeans_norm = True # Normalize the distances in each variable by the full range of values within each variable used for K-means clustering
make_cluster = False
compute_cluster_err = False
process_obj = False

err_threshold = [0.05, 0.10, 0.20, 0.30]

objnames = ['Energy','Water Quality','Flood']
n_obj = len(objnames)
altnames = ['Historical', 'WWTP', 'Basin', 'Leak']
n_alt = len(altnames)
parnames = ['HK_1','HK_2','HK_3','HK_4','HK_5','SS_1','SS_2','SS_3','SS_4','SS_5','SY_1','SY_2','SY_3','SY_4','SY_5','VANI_1','VANI_2','VANI_3','VANI_4','VANI_5','Q_1','Q_2','Q_3','LK_1','LK_2','LK_3','TWU_1','TWU_2','TWU_3','RCH_1','RCH_2','RCH_3','IN_1']
n_params = len(parnames)
parnameslong = ['Hydraulic Conductivity Lacustrine','Hydraulic Conductivity Alluvial','Hydraulic Conductivity Basaltic','Hydraulic Conductivity\nVolcaniclastic','Hydraulic Conductivity Andesitic','Specific Storage Lacustrine','Specific Storage Alluvial','Specific Storage Basaltic','Specific Storage\nVolcaniclastic','Specific Storage Andesitic','Specific Yield Lacustrine','Specific Yield Alluvial','Specific Yield Basaltic','Specific Yield\n Volcaniclastic','Specific Yield Andesitic','Vertical Anisotropy of\nHydraulic Conductivity Lacustrine','Vertical Anisotropy of\nHydraulic Conductivity Alluvial','Vertical Anisotropy of\nHydraulic Conductivity Basaltic','Vertical Anisotropy of\nHydraulic Conductivity Volcaniclastic','Vertical Anisotropy of\nHydraulic Conductivity Andesitic','Urban to Periurban\nPumping Multiplier 1990','Urban to Periurban\nPumping Multiplier 2000','Urban to Periurban\nPumping Multiplier 2010','Leak Multiplier 1990','Leak Multiplier 2000','Leak Multiplier 2010','Total Water Use\nMultiplier 1990','Total Water Use\nMultiplier 2000','Total Water Use\nMultiplier 2010','Recharge Multiplier\nUrban Land Use','Recharge Multiplier\nNatural Land Use','Recharge Multiplier\nWetland Land Use','Leak Infiltration Rate']
transform = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0]

## Observation Clustering
'''
This section calculates the sum of squared weighted residual for subsamples of the full observation set based on user inputs to a K-means clustering scheme that uses Latitude, Longitude, Observation hydraulic head, and time (optional)
'''
if my_rank==0:
    print('Creating K-means clusters', flush=True)
    # Set observation groups using k-means clustering based on observation characteristics (latitude, longitude, head elevation, and time [optional]). User chooses whether to use the absolute characteristic values (norm=False) or normalized based on the spread of the sample within each characteristic (norm=True)
    df, obsinfo, obsstats, obsformation = pltvm.process_hobs('20200403_100000', '00000', obsinfo_loaded=True)
    df.rename(columns={'absobserved': 'Z', 'time_series': 't'}, inplace=True)
    df['t_group'] = np.where(df['t']<=4018, 1, np.where((df['t']>4018) & (df['t']<=7671), 2, 3)) # Set model period grouping based on predetermined group limits
    df.set_index('obs_name', inplace=True)
    
    if make_cluster:
        group_label = sa.kmeans_obs(df, n_clusters=n_clusters, norm=kmeans_norm, time=kmeans_time)
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-cluster.csv')), group_label.astype(int), fmt='%i', delimiter=',')
        pltvm.plt_cluster(df, group_label, str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-cluster.png')))
        
    else:
        group_label = np.loadtxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-cluster.csv')), dtype='i4', delimiter=',')
    
    # Create combinations of clusters    
    lst = list(np.arange(0,n_clusters))
    combs = []
    
    for i in range(1, len(lst)+1):
        els = [list(x) for x in itertools.combinations(lst, i)]
        combs.extend(els)
    
    # Loop through all hobs output files to calculate Sum of Squared Weight Residual (soswr) based on the clusters defined in group_label
    hob_files = sorted(os.listdir(str(hobspath)))
    obj_files = sorted(os.listdir(str(objpath)))
else:
    obsstats = None
    group_label = None
    combs = None
    hob_files = None
    obj_files = None

obsstats = comm.bcast(obsstats, root=0)
group_label = comm.bcast(group_label, root=0)
combs = comm.bcast(combs, root=0)
n_combs = len(combs)
hob_files = comm.bcast(hob_files, root=0)
obj_files = comm.bcast(obj_files, root=0)

if compute_cluster_err:
    # Create a dictionary to store numpy array of error values for each k-means cluster
    kmc_err = {}
    for i in range(n_combs+1):
        kmc_err[i] = np.zeros(samples_per_proc)
    # Assign sample id for each error set
    kmc_err[0][:] = np.arange(0,samples_per_proc) + samples_per_proc*my_rank
    
    ## Error Calculation
    '''
    This section calculates the error using subsamples of the full observation set as determined by the clustering above. The error in this case study is determined as the sum of squared weighted residuals (soswr) with the weights determined from Mautner et al. 2020
    '''
    err_array = np.ones((samples_per_proc, n_clusters))*np.nan
    my_hob_files = hob_files[int(kmc_err[0][0]):int(kmc_err[0][kmc_err[0].shape[0]-1]+1)]
    print('Calculating error for samples ' + my_hob_files[0] + ' through ' + my_hob_files[len(my_hob_files)-1] + ' on ' + str(my_rank), flush=True)
    for s, s_file in enumerate(my_hob_files):   
        # Calculate soswr for sample based on full observation set
        hob_df = pd.read_fwf(hobspath.joinpath(s_file),widths=[22,19,22])
        heads_obs = [hob_df.columns.values.tolist()] + hob_df.values.tolist()
        
        # Calculate soswr for sample based on cluster observation set
        for c in range(n_clusters):
            heads = [heads_obs[i] for i, x in enumerate(np.append([True],[group_label==c])) if x]
            stats = obsstats.values[group_label==c]
            err_array[s,c], maxerror = mo.calculate_SOSWR(heads, stats) # assign soswr calculated for cluster to error array
        
    # Loop through all possible combinations of clusters and save combinations set
    combs_array = np.ones((n_combs,n_clusters))*np.nan
    for i, comb in enumerate(combs):
        for j in comb:
            kmc_err[i+1] += err_array[:,j]
        combs_array[i,:len(comb)] = np.array(comb)
        
    # Combine all error values across the processors
    comm.Barrier() # wait for all of them to finish
    kmc_err_all = comm.gather(kmc_err, root=0) # gather back to master node (0)

# the combined results only exist on master node
if my_rank==0:
    if compute_cluster_err:
        kmc_err_final = np.ones((n_samples, n_combs))*np.nan
        # Loop through list of dictionaries from each processor
        for d in kmc_err_all:
            # Loop through each column of error values for current processor samples
            for i in range(n_combs):
                kmc_err_final[d[0].astype(int),i] = d[i+1] # Use index entries from each dictionary to assign error values for each column
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-combs.csv')), combs_array, delimiter=',')
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-err.csv')), kmc_err_final, delimiter=',')
    
    else:
        kmc_err_final = np.loadtxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-err.csv')), delimiter=',')

## Management Objectives
'''
This section creates a dictionary for each management alternative that includes an n_samples by n_objectives matrix with the calculated value for each management objective
'''
if process_obj:
    # Create a dictionary to store numpy array of management objectives values for each sample
    sample_obj = {}
    for i in range(int(n_obj*n_alt+1)):
        sample_obj[i] = np.ones(samples_per_proc)*np.nan
    # Assign sample id for each error set
    sample_obj[0][:] = np.arange(0,samples_per_proc) + samples_per_proc*my_rank
    
    my_obj_files = obj_files[int(sample_obj[0][0]):int(sample_obj[0][sample_obj[0].shape[0]-1]+1)]
    print('Exctracting objectives for samples ' + my_obj_files[0] + ' through ' + my_obj_files[len(my_obj_files)-1] + ' on ' + str(my_rank), flush=True)
    for s, s_file in enumerate(my_obj_files):   
        # Exctract objective data for each sample
        tempobj = np.loadtxt(str(objpath.joinpath(s_file)), delimiter=',').flatten()
        
        # Loop through objectives for each scenario
        for o in range(tempobj.shape[0]):
            sample_obj[o+1][s] = tempobj[o]
    
    # Combine all obj values across the processors
    comm.Barrier() # wait for all of them to finish
    obj_all = comm.gather(sample_obj, root=0) # gather back to master node (0)

## Data Cleaning
'''
This section compiles error and objective data, removing any samples that contain nan in any objectives and preparing the data for sensitivity analysis
'''
# the combined results only exist on master node
if my_rank==0:
    if process_obj:
        print('Cleaning data', flush=True)
        obj_final = np.ones((n_samples, n_obj*n_alt))*np.nan
        for d in obj_all:
            for i in range(n_obj*n_alt):
                obj_final[d[0].astype(int),i] = d[i+1]
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-obj.csv')), obj_final, delimiter=',')
    else:
        obj_final = np.loadtxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-obj.csv')), delimiter=',')
        
    # Import parameters array
    params_array = np.loadtxt('params_' + folder + '.csv', delimiter=',')
#    c = 0
#    for column in params_array.T:
#        if transform[c] == 1:
#            params_array[:,c] = np.exp(column)
#        c += 1

    # Concatenate all datasets to be evaluated for Delta sensitivity, including parameter values
    sa_data_array = np.concatenate((np.reshape(np.arange(0,n_samples),(n_samples,1)), params_array, obj_final, kmc_err_final), axis=1)
    
    # Remove all samples that returned NaN for any objective values
    sa_data_array = sa_data_array[~np.isnan(sa_data_array).any(axis=1)]
    
    # Create dictionary of indices for datasets that are below the sample percent threshold for error of each cluster of observations
    sa_data = {}
    for i in err_threshold:
        thresh_sample = int(n_samples*i)
        sa_err_sort = np.zeros((thresh_sample, n_combs+1))
        sa_err_sort[:,n_combs] = i
        
        # Sort by the error for each cluster and record the indices of samples that are within the top i fraction of the total samples
        for j in range(n_combs):
            temp_sort = sa_data_array[:, 1+n_params+n_obj*n_alt+j].argsort()
            sa_err_sort[:,j] = temp_sort[:thresh_sample]
            
        
    np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-samples_evald.csv')), err_obj_array, delimiter=',') # Save filtered dataset
    
    # Generate problem set for parameters
    params = sa.gen_param_vals(err_obj_array.shape[0])
    params['values'] = err_obj_array[:,1:n_params+1]
    
else:
    df = None
    params = None

err_obj_array = comm.bcast(err_obj_array, root=0)
params = comm.bcast(params, root=0)

## Delta Sensitivity
'''
This section calculates the Delta sensitivity for each cluster soswr and each of the management objectives
'''
if my_rank < (n_obj*n_alt + n_combs):
    print('Calculating sensitivity for ' + str(my_rank), flush=True)
    delta_s = delta.analyze(params['problem'], params['values'], err_obj_array[:, my_rank + 1 + n_params])
    delta_s['column'] = my_rank
else:
    delta_s = None

# Combine all sensitivity values across the processors
comm.Barrier() # wait for all of them to finish
delta_all = comm.gather(delta_s, root=0) # gather back to master node (0)
# the combined results only exist on master node
if my_rank==0:
    print('Gathering sensitivity data', flush=True)
    delta_all_final = np.ones((n_params, n_obj*n_alt + n_combs))*np.nan
    for d in delta_all:
        if d is not None:
            delta_all_final[:, d['column']] = np.array(d['delta'])
        
    delta_obj_final = delta_all_final[:, :n_obj*n_alt]
    delta_clust_final = delta_all_final[:, n_obj*n_alt:]
    
    np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-deltas-obj.csv')), delta_obj_final, delimiter=',')
    np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-deltas-err.csv')), delta_clust_final, delimiter=',')
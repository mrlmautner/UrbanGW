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
import time

folder = '20200921_100000'
startpath = Path.cwd()
outputpath = startpath.joinpath('model_files').joinpath('output').joinpath('sa')
hobspath = outputpath.joinpath('hob').joinpath(folder)
objpath = outputpath.joinpath('obj').joinpath(folder)

comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()
num_proc = comm.Get_size()

n_clusters = 5
n_samples = 100000
samples_per_proc = int(n_samples / num_proc)

kmeans_time = False # Include time as a variable for K-means clustering
kmeans_norm = True # Normalize the distances in each variable by the full range of values within each variable used for K-means clustering
make_cluster = False # Boolean to create cluster if not already existing
compute_cluster_err = True # Process the error data for all clusters and model runs
process_obj = True # Process objective data for all clusters and model runs

err_threshold = [0.05, 0.10, 0.20, 0.30, 1] # Number of parameter sets for which to process error and objective data

objnames = ['Energy','Water Quality','Flood']
n_obj = len(objnames)
altnames = ['Historical', 'WWTP', 'Basin', 'Leak']
n_alt = len(altnames)
obj_alt = ['E-Hist','E-WWTP','E-Basin','E-Leak','W-Hist','W-WWTP','W-Basin','W-Leak','F-Hist','F-WWTP','F-Basin','F-Leak']
parnames = ['HK_1','HK_2','HK_3','HK_4','HK_5','SS_1','SS_2','SS_3','SS_4','SS_5','SY_1','SY_2','SY_3','SY_4','SY_5','VANI_1','VANI_2','VANI_3','VANI_4','VANI_5','Q_1','Q_2','Q_3','LK_1','LK_2','LK_3','TWU_1','TWU_2','TWU_3','RCH_1','RCH_2','RCH_3','IN_1']
n_params = len(parnames)
parnameslong = ['Hydraulic Conductivity Lacustrine','Hydraulic Conductivity Alluvial','Hydraulic Conductivity Basaltic','Hydraulic Conductivity\nVolcaniclastic','Hydraulic Conductivity Andesitic','Specific Storage Lacustrine','Specific Storage Alluvial','Specific Storage Basaltic','Specific Storage\nVolcaniclastic','Specific Storage Andesitic','Specific Yield Lacustrine','Specific Yield Alluvial','Specific Yield Basaltic','Specific Yield\n Volcaniclastic','Specific Yield Andesitic','Vertical Anisotropy of\nHydraulic Conductivity Lacustrine','Vertical Anisotropy of\nHydraulic Conductivity Alluvial','Vertical Anisotropy of\nHydraulic Conductivity Basaltic','Vertical Anisotropy of\nHydraulic Conductivity Volcaniclastic','Vertical Anisotropy of\nHydraulic Conductivity Andesitic','Urban to Periurban\nPumping Multiplier 1990','Urban to Periurban\nPumping Multiplier 2000','Urban to Periurban\nPumping Multiplier 2010','Leak Multiplier 1990','Leak Multiplier 2000','Leak Multiplier 2010','Total Water Use\nMultiplier 1990','Total Water Use\nMultiplier 2000','Total Water Use\nMultiplier 2010','Recharge Multiplier\nUrban Land Use','Recharge Multiplier\nNatural Land Use','Recharge Multiplier\nWetland Land Use','Leak Infiltration Rate']
transform = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0]

## Observation Clustering
'''
This section calculates the sum of squared weighted residual for subsamples of the full observation set based on user inputs to a K-means clustering scheme that uses Latitude, Longitude, Observation hydraulic head, and time (optional)
'''
if my_rank==0:
    # Set observation groups using k-means clustering based on observation characteristics (latitude, longitude, head elevation, and time [optional]). User chooses whether to use the absolute characteristic values (norm=False) or normalized based on the spread of the sample within each characteristic (norm=True)
    df, obsinfo, obsstats, obsformation = pltvm.process_hobs(folder, '00000', obsinfo_loaded=True)
    df.rename(columns={'absobserved': 'Z', 'time_series': 't'}, inplace=True)
    df['t_group'] = np.where(df['t']<=4018, 1, np.where((df['t']>4018) & (df['t']<=7671), 2, 3)) # Set model period grouping based on predetermined group limits
    df.set_index('obs_name', inplace=True)
    
    if make_cluster:
        print('Creating K-means clusters', flush=True)
        group_label = sa.kmeans_obs(df, n_clusters=n_clusters, norm=kmeans_norm, time=kmeans_time)
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-cluster.csv')), group_label.astype(int), fmt='%i', delimiter=',')
        pltvm.plt_cluster(df, group_label, str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-cluster.png')))
        
    else:
        print('Loading K-means clusters', flush=True)
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
    kmc_err_all = comm.gather(kmc_err, root=0) # gather back to conductor node (0)

# the combined results only exist on conductor node
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
    obj_all = comm.gather(sample_obj, root=0) # gather back to conductor node (0)

## Data Cleaning
'''
This section compiles error and objective data, removing any samples that contain nan in any objectives and preparing the data for sensitivity analysis
'''
# the combined results only exist on conductor node
if my_rank==0:
    if process_obj:
        print('Cleaning data', flush=True)
        obj_final = np.ones((n_samples, n_obj*n_alt))*np.nan
        for d in obj_all:
            for i in range(n_obj*n_alt):
                obj_final[d[0].astype(int),i] = d[i+1]
        np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + '-obj.csv')), obj_final, delimiter=',')
    else:
        obj_final = np.loadtxt(str(outputpath.joinpath('analysis').joinpath(folder + '-obj.csv')), delimiter=',')
        
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
    sa_data_sort = {}
    sa_data_index = np.ones((len(err_threshold)*n_combs*n_obj*n_alt,3))#np.ones((len(err_threshold)*n_obj*n_alt,3))#
    a = 0 # index for error threshold
    b = 0 # index for comb
    for i in err_threshold:
        thresh_sample = int(n_samples*i)
        sa_err_sort = np.zeros((thresh_sample, n_combs+1))
        sa_err_sort[:,n_combs] = i
        sa_data_index[a:a+n_combs*n_obj*n_alt, 0] = i # index that represents error threshold
#        sa_data_index[a:a+n_obj*n_alt, 0] = i
        a += n_combs*n_obj*n_alt
#        a += n_obj*n_alt
        
        # Sort by the error for each cluster and record the indices of samples that are within the top i fraction of the total samples
        for j in range(n_combs):
            temp_sort = sa_data_array[:, 1+n_params+n_obj*n_alt+j].argsort()
            sa_err_sort[:,j] = temp_sort[:thresh_sample]
            sa_data_index[b:b+n_obj*n_alt, 1] = j # index that represents cluster combination
            sa_data_index[b:b+n_obj*n_alt, 2] = np.arange(n_obj*n_alt) # index that represents objective/error function
            b += n_obj*n_alt
#            if j == n_combs-1:
#                sa_data_index[b:b+n_obj*n_alt, 1] = j
#                sa_data_index[b:b+n_obj*n_alt, 2] = np.arange(n_obj*n_alt)
#                b += n_obj*n_alt
            
        sa_data_sort[i] = sa_err_sort
    
    np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-samples_evald.csv')), sa_data_array, delimiter=',') # Save filtered dataset
    np.savetxt(str(outputpath.joinpath('analysis').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters) + '-threshold_samples.csv')), np.concatenate(tuple(list(sa_data_sort.values())), axis=0), delimiter=',') # Save filtered dataset
    
else:
    sa_data_sort = None
    sa_data_index = None
    sa_data_array = None

sa_data_sort = comm.bcast(sa_data_sort, root=0)
sa_data_index = comm.bcast(sa_data_index, root=0)
sa_data_array = comm.bcast(sa_data_array, root=0)

## Delta Sensitivity
'''
This section calculates the Delta (delta) and Sobol 1st order (S1) sensitivity for each cluster soswr and each of the management objectives
'''
# Define conductor/player routine to calculate delta sensitivity
WORKTAG = 1
DIETAG = 0

def conductor(comm):
    status = MPI.Status()
    
    sarun = 0
    delta_data = {}
    S1_data = {}
    for i in err_threshold:
        delta_data[i] = np.zeros((n_combs, n_params, n_obj*n_alt))
        S1_data[i] = np.zeros((n_combs, n_params, n_obj*n_alt))
    
    # Seed the players, send one unit of work to each player (rank)
    for rank in range(1, int(min(num_proc, sa_data_index.shape[0]))):
        comm.send(sarun, dest=rank, tag=WORKTAG)
        sarun += 1
    
    # Loop over getting new work requests until there is no more work to be done
    while True:
        
        # Receive results from a player
        player_data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        player_index = player_data[0]
        delta_vals = np.array(player_data[1]['delta'])
        S1_vals = np.array(player_data[1]['S1'])
        delta_data[player_index[0]][int(player_index[1]),:,int(player_index[2])] = delta_vals
        S1_data[player_index[0]][int(player_index[1]),:,int(player_index[2])] = S1_vals
        
        if sarun >= sa_data_index.shape[0]:
            break

        # Send the player a new work unit
        comm.send(sarun, dest=status.Get_source(), tag=WORKTAG)
        sarun += 1
            
    # No more work to be done, receive all outstanding results from players
    print('Receiving outstanding work', flush=True)
    for rank in range(1, int(min(num_proc-1, sa_data_index.shape[0]))): 
        player_data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        player_index = player_data[0]
        delta_vals = np.array(player_data[1]['delta'])
        S1_vals = np.array(player_data[1]['S1'])
        delta_data[player_index[0]][int(player_index[1]),:,int(player_index[2])] = delta_vals
        S1_data[player_index[0]][int(player_index[1]),:,int(player_index[2])] = S1_vals

    # Tell all the players to exit by sending an empty message with DIETAG
    print('Killing processers', flush=True)
    for rank in range(1, num_proc):
        comm.send(0, dest=rank, tag=DIETAG)

    # Save the sum-of-squared, weighted residual error
    print('Saving sensitivity values', flush=True)
    sensitivity_loc = outputpath.joinpath('sens').joinpath(folder + ('-norm' if kmeans_norm else '') + ('-t' if kmeans_time else '') + '-' + str(n_clusters))
    for i in delta_data.keys():
        for j in range(n_combs):
            try:
                np.savetxt(str(sensitivity_loc.joinpath('delta-' + str(i) + '-' + '{:02d}'.format(j) + '.csv')), delta_data[i][j,:,:], delimiter=',')
                np.savetxt(str(sensitivity_loc.joinpath('S1-' + str(i) + '-' + '{:02d}'.format(j) + '.csv')), S1_data[i][j,:,:], delimiter=',')
            except:
                sensitivity_loc.mkdir(exist_ok=True)
                np.savetxt(str(sensitivity_loc.joinpath('delta-' + str(i) + '-' + '{:02d}'.format(j) + '.csv')), delta_data[i][j,:,:], delimiter=',')
                np.savetxt(str(sensitivity_loc.joinpath('S1-' + str(i) + '-' + '{:02d}'.format(j) + '.csv')), S1_data[i][j,:,:], delimiter=',')

def player(comm):
    my_rank = comm.Get_rank()
    status = MPI.Status()

    while True:
        # Receive a message from the conductor
        sarun = int(comm.recv(source=0, tag=MPI.ANY_TAG, status=status))

        # Check the tag of the received message
        if status.Get_tag() == DIETAG: break 

        timerun = time.time()
        print('Calculating sensitivity for model '+ str(sarun) + ' on ' + str(my_rank), flush=True)
        
        sarun_index = sa_data_index[sarun,:]
        
        # Generate problem set for parameters
        params = sa.gen_param_vals(sa_data_sort[sarun_index[0]].shape[0])
        params['values'] = sa_data_array[sa_data_sort[sarun_index[0]][:,int(sarun_index[1])].astype('int'),1:n_params+1]
        
        # Calculate sensitivity of the column selected using the threshold, cluster, and objective/error selection indicated by sa_data_index[0], sa_data_index[1], and sa_data_index[2] respectively
        try:
            sens_dict = delta.analyze(params['problem'], params['values'], sa_data_array[sa_data_sort[sarun_index[0]][:,int(sarun_index[1])].astype('int'), 1+n_params+int(sarun_index[2])])
        except:
            sens_dict = np.ones((n_obj*n_alt))*np.nan
            
        print('Finished model number '+ str(sarun) + ' on ' + str(my_rank) + ' in ' + str(time.time() - timerun) + ' seconds', flush=True)

        # Send the result back
        comm.send([sarun_index, sens_dict], dest=0, tag=0)
        
if my_rank == 0:
    conductor(comm)
else:
    player(comm)
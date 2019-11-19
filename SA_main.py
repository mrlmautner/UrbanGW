import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
import numpy as np
from mpi4py import MPI
import time

tot_samples = 10080
soswrlim = 1000000
safolder = '20191028_10000'

# MODFLOW File location
exefile = '/zeolite/mmautner/MODFLOW/MF2005.1_12u/make/mf2005'

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_proc = comm.Get_size()

if rank == 0:
    params = sa.gen_param_vals(tot_samples)
    np.savetxt('params_' + safolder + '.csv', params['values'], delimiter=',')
else:
    params = None
params = comm.bcast(params, root=0)

if rank == 0:
    # Load alternatives
    altpath = Path.cwd() / 'model_files' / 'optimization_data' / 'decisions' / 'alternatives.csv'
    alternatives = {}
    altkeys = ['names','wwtps','basins','leakrepair','labels']
    with open(altpath) as a:
        lines = a.read().splitlines() 
        for i, line in enumerate(lines):
            templist = line.split(',')
            alternatives[altkeys[i]] = templist
else:
    alternatives = None
alternatives = comm.bcast(alternatives, root=0)

samples_per_proc = tot_samples/num_proc
for i in range(int(samples_per_proc)):
    sarun = int(rank * samples_per_proc + i)
    timerun = time.time()
    print('Running model number '+ str(sarun) + ' on ' + str(rank))
    
    try:
        objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=False)
    except:
        objectives = np.nan
        
    print('Finished model number '+ '{:05d}'.format(sarun) + ' on ' + str(rank) + ' in ' + str(time.time() - timerun) + ' seconds')
    
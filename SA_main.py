import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
from mpi4py import MPI

tot_samples = 10
soswrlim = 10000000
safolder = '20191022'

# MODFLOW File location
exefile = Path.cwd() / 'gwscripts' / 'gwmodel' / 'make' / 'mf2005'

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_proc = comm.Get_size()

if rank == 0:
    params = sa.gen_param_vals(tot_samples)
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
for i in range(samples_per_proc):
    sarun = rank * samples_per_proc + i

    objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=False)
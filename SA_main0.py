import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
#from mpi4py import MPI

# MODFLOW File location
exefile = Path.cwd() / 'gwscripts' / 'gwmodel' / 'make' / 'mf2005'

# Load alternatives
altpath = Path.cwd() / 'model_files' / 'optimization_data' / 'decisions' / 'alternatives.csv'
alternatives = {}
altkeys = ['names','wwtps','basins','leakrepair','labels']
with open(altpath) as a:
    lines = a.read().splitlines() 
    for i, line in enumerate(lines):
        templist = line.split(',')
        alternatives[altkeys[i]] = templist

nsamples = 1000
soswrlim = 1000000000
params = sa.gen_param_vals(nsamples)
safolder = '20191017'
sarun = 0

objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=True)
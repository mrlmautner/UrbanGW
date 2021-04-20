import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
import numpy as np
from mpi4py import MPI
import time

tot_samples = 100000
soswrlim = 2000000
safolder = '20200921_100000'
numSets = 20
currentSet = 7
sarun = int((tot_samples/numSets)*currentSet)

# vvvvvvvvvvvv Do not touch below this line vvvvvvvvvvvv
# MODFLOW File location
exefile = '/zeolite/mmautner/MODFLOW/MF2005.1_12u/make/mf2005'
comm = MPI.COMM_WORLD
my_rank = comm.Get_rank()

if my_rank == 0:
    # Create or load params file
    paramspath = 'params_' + safolder + '.csv'
    try:
        with open(paramspath) as a:
            lines = a.read().splitlines() 
            templist = []
            for i, line in enumerate(lines):
                templist.append([float(x) for x in line.split(',')])
        params = sa.gen_param_vals(tot_samples)
        params['values'] = np.array(templist)
    except:
        params = sa.gen_param_vals(tot_samples)
        np.savetxt(paramspath, params['values'], delimiter=',')
else:
    params = None
params = comm.bcast(params, root=0)

if my_rank == 0:
    # Load alternatives
    altpath = 'model_files/optimization_data/decisions/alternatives.csv'
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
    
WORKTAG = 1
DIETAG = 0

def master(comm):
    num_procs = comm.Get_size()
    status = MPI.Status()
    
    sarun = int((tot_samples/numSets)*currentSet)

    # Seed the slaves, send one unit of work to each slave (rank)
    for rank in range(1, int(min(num_procs,1+tot_samples/numSets))):
        comm.send(sarun, dest=rank, tag=WORKTAG)
        sarun += 1
        
    # Loop over getting new work requests until there is no more work to be done
    while True:
            
        # Receive results from a slave
        objectives = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

        if sarun >= int((tot_samples/numSets)*(currentSet+1)):
            break
    
        # Send the slave a new work unit
        comm.send(sarun, dest=status.Get_source(), tag=WORKTAG)
        sarun += 1

    # Tell all the slaves to exit by sending an empty message with DIETAG
    for rank in range(1, num_procs):
        comm.send(0, dest=rank, tag=DIETAG)
  
def slave(comm):
    status = MPI.Status()
    
    while True:
        # Receive a message from the master
        sarun = int(comm.recv(source=0, tag=MPI.ANY_TAG, status=status))
    
        # Check the tag of the received message
        if status.Get_tag() == DIETAG:
            break 
    
        # Do the work
        timerun = time.time()
        print('Running model number '+ str(sarun) + ' on ' + str(my_rank), flush=True)
            
        #objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=False)
        try:
            objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=False)
        except:
            objectives = np.nan
        
        print('Finished model number '+ '{:05d}'.format(sarun) + ' on ' + str(my_rank) + ' in ' + str(time.time() - timerun) + ' seconds', flush=True)
    
        # Send the result back
        comm.send(objectives, dest=0, tag=0)
  
def main():
    if my_rank == 0:
        master(comm)
    else:
        slave(comm)

if __name__ == '__main__':
    main()
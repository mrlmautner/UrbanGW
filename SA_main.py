import model_functions as mf
import gwscripts.sensitivityanalysis.satools as sa
from pathlib import Path
import numpy as np
from mpi4py import MPI
import time

tot_samples = 400
soswrlim = 1000000
safolder = '20191230_400'

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

WORKTAG = 1
DIETAG = 0

def master(comm):
    num_procs = comm.Get_size()
    status = MPI.Status()

    sarun = 0
    # Seed the slaves, send one unit of work to each slave (rank)
    for rank in range(1, num_procs):
        comm.send(sarun, dest=rank, tag=WORKTAG)
        sarun += 1
    
    # Loop over getting new work requests until there is no more work to be done
    while True:
        
        # Receive results from a slave
        objectives = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)

        # Send the slave a new work unit
        comm.send(sarun, dest=status.Get_source(), tag=WORKTAG)
        sarun += 1
        if sarun == tot_samples:
            break
    
    # No more work to be done, receive all outstanding results from slaves
    for rank in range(1, num_procs): 
        objectives = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        #process_result(result)

    # Tell all the slaves to exit by sending an empty message with DIETAG
    for rank in range(1, num_procs):
        comm.send(0, dest=rank, tag=DIETAG)


def slave(comm):
    my_rank = comm.Get_rank()
    status = MPI.Status()

    while True:
        # Receive a message from the master
        sarun = int(comm.recv(source=0, tag=MPI.ANY_TAG, status=status))

        # Check the tag of the received message
        if status.Get_tag() == DIETAG: break 

        # Do the work
        timerun = time.time()
        print('Running model number '+ str(sarun) + ' on ' + str(my_rank))
        
        try:
            objectives = mf.SA_mode(alternatives=alternatives, params=params, exefile=exefile, safolder=safolder, sarun=sarun, soswrlim=soswrlim, verbose=False)
        except:
            objectives = np.nan
            
        print('Finished model number '+ '{:05d}'.format(sarun) + ' on ' + str(my_rank) + ' in ' + str(time.time() - timerun) + ' seconds')

        # Send the result back
        comm.send(objectives, dest=0, tag=0)
        

def main():
    comm = MPI.COMM_WORLD
    my_rank = comm.Get_rank()
    my_name = MPI.Get_processor_name()
    #comm.Barrier()
    #start = MPI.Wtime()
    
    if my_rank == 0:
        master(comm)
    else:
        slave(comm)
    
    #comm.Barrier()
    #end = MPI.Wtime()
    #print 'time:', end - start

if __name__ == '__main__':
    main()
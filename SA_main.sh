#!/bin/bash -l
#SBATCH -n 168                  # Total number of processors to request
#SBATCH --distribution=plane=8  # Distribute in blocks of 8 to each node
#SBATCH -p aqua                 # Partition option of AQUA is aqua
#SBATCH -t 200:00:00            # Run time (hh:mm:ss) - 24 hours
#SBATCH -M mmautner@ucdavis.edu              # address for email notification
#SBATCH --mail-type=ALL         # email at Begin and End of job

# Add my directory to the path to point to python.
export PATH=/zeolite/mmautner/miniconda3/bin:$PATH

# To point to the correct MPI in my local folder
module purge 

mpirun --bind-to none python SA_main.py
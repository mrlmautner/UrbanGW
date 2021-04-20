#!/bin/bash -l
#SBATCH -n 3                    # Total number of processors to request
#SBATCH -p aqua                 # Partition option of AQUA is aqua
#SBATCH -t 200:00:00            # Run time (hh:mm:ss) - 24 hours
#SBATCH --mail-user=mmautner@ucdavis.edu              # address for email notification
#SBATCH --mail-type=ALL         # email at Begin and End of job

# Add my directory to the path to point to python.
export PATH=/zeolite/mmautner/miniconda3/bin:$PATH

# To point to the correct MPI in my local folder
module purge 

mpirun --bind-to none python SA_main_missing.py
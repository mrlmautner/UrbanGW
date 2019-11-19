#!/bin/bash -l
#SBATCH -n 1            # Total number of processors to request (32 cores per node)
#SBATCH -p aqua
#SBATCH -t 200:00:00        # Run time (hh:mm:ss) - 24 hours
#SBATCH --mail-user=mmautner@ucdavis.edu              # address for email notification
#SBATCH --mail-type=ALL                  # email at Begin and End of job

# you can define your own variables. Access them later with dollar signs ($DIR, etc.)

# IMPORTANT: Python3 is installed in my directory. Add it to the path.
export PATH=/zeolite/mmautner/miniconda3/bin:$PATH

module purge

python SA_main0.py
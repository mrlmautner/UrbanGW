#!/bin/bash -l
#SBATCH -n 200                  # Total number of processors to request
#SBATCH -p high                # Queue name high/med/low
#SBATCH -t 200:00:00           # Run time (hh:mm:ss) - 24 hours
#SBATCH --mail-user=mmautner@ucdavis.edu              # address for email notification
#SBATCH --mail-type=ALL        # email at Begin and End of job

# you can define your own variables. Access them later with dollar signs ($DIR, etc.)
GDIR=/group/hermangrp

# IMPORTANT: Python3/Pyomo/CBC solver are all installed in group directory. Add it to the path.
export PATH=$GDIR/miniconda3/bin:$GDIR/cbc/bin:/home/mmautner/MODFLOW/MF2005.1_12u/make:$PATH

mpirun python process_SA.py
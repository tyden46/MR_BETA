#!/bin/bash

#SBATCH --time=00:30:00   # walltime
#SBATCH --ntasks=16   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
module load python/2.7.16
python BEAR_release/scripts/generate_reads2.py -r merged.fna -a coreMedianAbundance.txt -o 10Million-2.fna -t 10000000 -l 200

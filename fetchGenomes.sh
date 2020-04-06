#!/bin/bash

#SBATCH --time=00:15:00   # walltime
#SBATCH --ntasks=4   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
SAVEIFS=$IFS   # Save current IFS
IFS=$'\n'      # Change IFS to new line
fileName=$1
mapfile -t list < $fileName
declare -i counter=1
for i in "${list[@]}"
do
  ncbi-acc-download --format fasta --out temp-myGenome${counter}.fna $i
  counter=$(expr $counter + 1)
done
cat temp* > myGenomes.fna
rm temp-*
IFS=$SAVEIFS   # Restore IFS

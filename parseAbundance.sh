#!/bin/bash

#SBATCH --time=14:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "parseAbundance"       # job name
SAVEIFS=$IFS   # Save current IFS
IFS=$'\n'      # Change IFS to new line
fileName=$1
fileLength=$(wc -l $fileName)
numLines=${fileLength%% *}
numReads=$(expr $numLines / 2)
grep -E -o '^[>]{1}[^[:space:]]+' $fileName | sort --unique > tempList.txt
mapfile -t list < tempList.txt
rm tempList.txt
for index in {0..71}
do
  id="${list[index]}"
  id="${id:1}"
  numForGenome=$(grep -c "${list[index]}" $fileName)
  proportion=$(bc <<< "scale=8; $numForGenome / $numReads")
  printf "%s\t%s\n" $id $proportion
done
IFS=$SAVEIFS   # Restore IFS

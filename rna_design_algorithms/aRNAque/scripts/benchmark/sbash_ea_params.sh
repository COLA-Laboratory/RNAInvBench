#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=65536
#SBATCH --partition=bdw
#SBATCH --mail-user=nonosaha@mis.mpg.de 
#SBATCH --mail-type=ALL

./mut_benchmark.sh $1 $2 

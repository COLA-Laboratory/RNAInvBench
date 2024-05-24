#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=65536
#SBATCH --partition=bdw
#SBATCH --mail-user=nonosaha@mis.mpg.de 
#SBATCH --mail-type=ALL
#SBATCH --time-limit=ULIMITED
python ea_one2.py -n 10 -g 10000  -t $1 -c 1.6 --job=20 
 

#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=65536
#SBATCH --partition=bdw
#SBATCH --mail-user=nonosaha@mis.mpg.de 
#SBATCH --mail-type=ALL

len=`expr length $1`;
step=`bc <<< "scale=4; 0.5/$len"`
for mu in `seq $step $step 0.2 `;	
	do 
		`python aRNAque_pk.py -t $1 -ft "ip" -n 100 -g 200 -msf 1 -mu $mu --job 20 -sm "F" >>../data/pkbase/pk_tuning/pk__mu"$len".out`;
	      	echo "solved for mu=: $mu target = $1"; 
	done


#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=65536
#SBATCH --partition=ivy
#SBATCH --mail-user=nonosaha@mis.mpg.de 
#SBATCH --mail-type=ALL

len=`expr length $1`;
step=0.1
for c in `seq 1 $step 2 `;	
	do 
		`python aRNAque_pk.py -t $1 -ft "ip" -n 100 -g 200 -msf 1 -c $c --job 20 -sm "F" >>../data/pkbase/pk_tuning/pk__c"$len".out`;
	      	echo "solved for c=: $c target = $1"; 
	done


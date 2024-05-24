#!/bin/bash 
for bp in `cat best_params`; 
	do 
		`python clean_aRNAque.py -t "$1"  -n 100 -g $2 -msf 1  --job 20 -bp="$bp" -c 1.6 --hairpin_boosting >> log_mut_levy_bench2.txt`;
	      	echo "solved for target: $1 and mut param : $bp";
	      	
	
	done


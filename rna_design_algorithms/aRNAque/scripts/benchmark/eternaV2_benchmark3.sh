#!/bin/bash 
for ss in `cat ../data/eterna/eterna2_unsolved_sorted.csv| cut -f 2 -d ","`; 
	do 
		`python aRNAque_new.py -t $ss  -n 100 -g 10000 -msf 1  --job 20 -c 1.6 -bp="GC1" >>../data/eterna/EternaV2_LevyGC1_10000.out`;
	      	echo "solved for target: $ss"; 
	
	done


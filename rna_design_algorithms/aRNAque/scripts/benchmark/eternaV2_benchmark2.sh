#!/bin/bash 
for ss in `cat ../data/eterna/unsolved_arnaqueGC1.txt| cut -f 2 -d ","`; 
	do 
		`python aRNAque_new.py -t $ss  -n 100 -g 200 -msf 1  --job 20 -c 1.6 -bp="GC" >>../data/eterna/EternaV2_LevyGC.out`;
	      	echo "solved for target: $ss"; 
	
	done


#!/bin/bash
for ss in `cat ../../data/Eterna100/V2/Eterna100.txt| cut -f 2 -d ","`;
	do
		`python ../../src/aRNAque.py -t $ss  -n 100 -g 10000 -msf 1  --job 1 -bp="GC4" >>../../data/eterna/EternaV2_OP4.out`;
	      	echo "solved for target: $ss";

	done

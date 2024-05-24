#!/bin/bash
for ss in `cat ../data/pkbase/rest.csv| cut -f 4 -d ","`;
	do
		`python aRNAque.py -t $ss -ft 'hk' -n 100 -g 200 -msf 1 -sm 'F' --job 20 -c 1.6 >>../data/pk_hotknotsLevy.out`;
	      	echo "solved for target: $ss";

	done

#!/bin/bash 

for target in `cat full_sorted_pk.csv| cut -f 3 -d ','`; 
	do 
		echo $target, $1 >> antaRNA_GC$1.out;
	       	`time python antaRNA_HKwrapper.py -t $target -nr 20 -it 5000 -tGC $1 --ft "IPKnot" >> antaRNA_GC$1_cputime.out` ;
      	done



#!/bin/bash 

for target in `cat full_sorted_pk$1.csv| cut -f 3 -d ','`; 
	do 
		echo $target, $1 >> antaRNA_GCHotknot$1.out;
	       	`time python antaRNA_HKwrapper.py -t $target -nr 20 -it 5000 -tGC .0 --ft "pKiss" >> antaRNA_GC_Hotknot_CPU$1.out` ;
      	done



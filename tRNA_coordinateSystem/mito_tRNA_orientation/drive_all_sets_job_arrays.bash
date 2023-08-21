#!/bin/bash

SOURCE=/scratch/hassan.as/mito.discovery.2021b

cd $SOURCE
for paramset in `ls S-0.3.Lpr-0.10.Cpr-0.2.Lpo-0.10.Cpo-0.8.t-1.0.r-APPE.tm-AAPP.Ab.expanded.fixed.L1-tRNA-8 -d`
do
        cd $SOURCE/$paramset
        for Set in `ls set* -d`
	do
		cd $SOURCE/$paramset		
		for run in `find $Set -maxdepth 1 -mindepth 1 -type d -printf '%f\n'`
        	do
                	cd $SOURCE/$paramset/$Set/$run
                	if [ -e STOP ]
               		then
                        	mkdir -p $SOURCE/Principal_axes_analysis/$paramset/$Set/$run
                        	for xtc in run*.xtc
				do
					echo "$SOURCE/$paramset/$Set/$run/$xtc" >> $SOURCE/Principal_axes_analysis/$paramset/$Set/$run/XTClist
				done
				cd $SOURCE/Principal_axes_analysis/$paramset/$Set/$run
				len=`wc -l XTClist | awk '{print $1}'`  #length of the array of xtc files
				let l=$len-1
				
				export paramset Set run 
				sbatch --array=0-$l $SOURCE/principal_axes.bash
			fi
			
        	done
	done
done


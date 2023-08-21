#!/bin/bash

cd Principal_axes_analysis
for i in S-0.3.Lpr-0.10.Cpr-0.2.Lpo-0.10.Cpo-0.8.t-1.0.r-APPE.tm-AAPP.Ab.expanded.fixed.L1-tRNA-8
do	

	cd $i
        for k in `ls`
        do
                cd $k
		for j in `ls`
	        do
                	cd $j
                	rm slurm*
			rm XTClist
			
			#tRNA principal axes processing:
			rm *paxis2* *paxis3* *moi*        
			cat `ls AtRNA_acceptor_paxis1_part*` > AtRNA_acceptor_paxis1.xvg
			cat `ls PtRNA_acceptor_paxis1_part*` > PtRNA_acceptor_paxis1.xvg
			cat `ls AtRNA_ASL_paxis1_part*` > AtRNA_ASL_paxis1.xvg
			cat `ls PtRNA_ASL_paxis1_part*` > PtRNA_ASL_paxis1.xvg
			
			rm *_paxis1_part*
			
#			/scratch/hassan.as/mito.discovery.2021b/calculate_principalAxes_angles.bash AtRNA_acceptor PtRNA_acceptor
#			/scratch/hassan.as/mito.discovery.2021b/calculate_principalAxes_angles.bash AtRNA_ASL PtRNA_ASL

#			rm *paxis*
#			rm *_X_*
			
#			paste AtRNA_acceptor_PtRNA_acceptor_angles AtRNA_ASL_PtRNA_ASL_angles | awk '{print $1,$2,$3,$5,$6}' > tRNA_axes_angles.xvg
#               	rm AtRNA_acceptor_PtRNA_acceptor_angles AtRNA_ASL_PtRNA_ASL_angles

#                	touch README
#                	echo -e "tRNA_axes_angles are: tRNA_acceptors_plane_vector-with-REF tRNA_acceptor_angle tRNA_ASL_plane_vector-with-REF tRNA_ASL_angle " > README
        		cd ..
        	done
		cd ..
        done
        cd ..
	echo "DONE: $i"
done


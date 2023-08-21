#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH -p redwood
#SBATCH --cpus-per-task=1
##SBATCH --exclusive
#SBATCH --job-name=tRNA_Paxis
#SBATCH --mem=40Gb
#SBATCH --output=slurm.out-%A-%a 
#SBATCH --error=slurm.err-%A-%a

module load gcc/6.4.0 openmpi/3.1.2
GMX=/home/whitford/BIN/nd/gmx-5.1.4_gnu6_openmpi3.1.2_redwood/bin/gmx_mpi

#XTC=$1
xtclist=($(cat XTClist))
XTC=${xtclist[${SLURM_ARRAY_TASK_ID}]}
echo $XTC

number=`echo $XTC | sed 's/^.*run.part00//' | sed 's/.xtc//'`  #extract run file number
#NDX=/scratch/hassan.as/mito.discovery.2021b/index_for_principalAxes_withL1.ndx
#GRO=/scratch/hassan.as/mito.discovery.2021b/minimized/S-0.3.Lpr-0.10.Cpr-0.2.Lpo-0.10.Cpo-0.8.t-1.0.r-APPE.tm-AAPP.Ab.expanded.fixed.L1-tRNA-8.min.gro


#echo 0 | $GMX principal -f $XTC -s $GRO -n $NDX -a1 AtRNA_ASL_paxis1_part_$number -a2 AtRNA_ASL_paxis2_part_$number -a3 AtRNA_ASL_paxis3_part_$number -om AtRNA_ASL_moi_part_$number -xvg none 

#echo 1 | $GMX principal -f $XTC -s $GRO -n $NDX -a1 AtRNA_acceptor_paxis1_part_$number -a2 AtRNA_acceptor_paxis2_part_$number -a3 AtRNA_acceptor_paxis3_part_$number -om AtRNA_acceptor_moi_part_$number -xvg none 

#echo 2 | $GMX principal -f $XTC -s $GRO -n $NDX -a1 PtRNA_ASL_paxis1_part_$number -a2 PtRNA_ASL_paxis2_part_$number -a3 PtRNA_ASL_paxis3_part_$number -om PtRNA_ASL_moi_part_$number -xvg none

#echo 3 | $GMX principal -f $XTC -s $GRO -n $NDX -a1 PtRNA_acceptor_paxis1_part_$number -a2 PtRNA_acceptor_paxis2_part_$number -a3 PtRNA_acceptor_paxis3_part_$number -om PtRNA_acceptor_moi_part_$number -xvg none 


################# Rotation angles:
SOURCE=/scratch/hassan.as/mito.discovery.2021b
path=/scratch/hassan.as/mito.discovery.2021b/Principal_axes_analysis/$paramset/$Set/$run   #exported variables from drive_all_sets_job_array.bash script

sed "s|FILE|$XTC|g; s/NUM/$number/g" $SOURCE/script_rot.vmd > $path/tmp.script_rot_$number\.vmd   #in the first sed command, I used the separator | to deal with the / that is present in $XTC
vmd -dispdev text -e $path/tmp.script_rot_$number\.vmd
rm $path/tmp.script_rot_$number\.vmd


#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH -p redwood
#SBATCH --cpus-per-task=1
##SBATCH --exclusive
#SBATCH --job-name=tRNA_rot
#SBATCH --mem=10Gb
#SBATCH --output=slurm.out-%A-%a 
#SBATCH --error=slurm.err-%A-%a

SOURCE=/scratch/hassan.as/mito_tRNA_orientation
vmd -dispdev text -e $SOURCE/tRNA_rotations.vmd

awk '{print NR,$1,$2}' AtRNA_orientation_traj > atrna.xvg
awk '{print NR,$1,$2}' PtRNA_orientation_traj > ptrna.xvg

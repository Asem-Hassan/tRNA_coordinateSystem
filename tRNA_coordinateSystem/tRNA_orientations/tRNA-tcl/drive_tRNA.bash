#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH -p short
#SBATCH --cpus-per-task=1
##SBATCH --exclusive
#SBATCH --job-name=tRNA_orien
#SBATCH --mem=10Gb
#SBATCH --output=slurm.out-%A-%a 
#SBATCH --error=slurm.err-%A-%a

SOURCE=/scratch/hassan.as/tRNA_orientations/tRNA-tcl
Matrices_dir=/scratch/hassan.as/tRNA_orientations/OUTPUT_Matrices

tRNApdbFile=$1
tRNAchainIDfinal=$2
NAME=$3
matrix_file=$4

MLARGE=`grep "Mlarge:" $Matrices_dir/$matrix_file | sed 's/Mlarge://'`
MBODY=`grep "Mbody:" $Matrices_dir/$matrix_file | sed 's/Mbody://'`

echo $MLARGE 
echo $MBODY

sed "s|PDBFILE|$tRNApdbFile|g; s/tRNAchainID/$tRNAchainIDfinal/g; s/MLARGE/$MLARGE/g; s/MBODY/$MBODY/g; s/NAME/$NAME/g" $SOURCE/tRNA.vmd > $SOURCE/OUTPUT/tRNA_tmp_$NAME.vmd

vmd -dispdev text -e $SOURCE/OUTPUT/tRNA_tmp_$NAME.vmd
rm $SOURCE/OUTPUT/tRNA_tmp_$NAME.vmd

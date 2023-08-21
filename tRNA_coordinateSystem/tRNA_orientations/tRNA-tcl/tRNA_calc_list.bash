#!/bin/bash

list=$1

SOURCE=/scratch/hassan.as/tRNA_orientations
FINALfile=$SOURCE/test_tRNA_list/tRNA_all_gt50
PDB_dir=/scratch/hassan.as/tRNA_orientations/PDB
Matrices_dir=/scratch/hassan.as/tRNA_orientations/OUTPUT_Matrices

#mkdir OUTPUT

pdblist=$( cat $list )

for pdb in $pdblist
do
	let tRNA_number=0
	numoflines=`cat $FINALfile | grep $pdb | wc -l`

	for line_nm in `seq 1 $numoflines`
	do
		line=`cat $FINALfile | grep $pdb | sed -n $line_nm\p`                #print line number $line_nm
		echo $line
		numoffields=`echo "$line" | awk 'BEGIN{FS=","}{for(i=6;i<=10;++i){if($i==""){print i-1;exit}}}'`
#		let numoftRNAs=$numoffields-5
		for tRNA in `seq 6 $numoffields`    
		do	
			tRNAchainID=`echo "$line" | awk -v col=$tRNA 'BEGIN{FS=","}{print $col}' | sed 's/ //g'`
			echo $tRNAchainID	
			pdb_lower=`echo "$pdb" | awk '{print tolower($0)}'`		
			tRNAchainIDfinal=`awk -v tRNA=$tRNAchainID '{if($2==""){pdbfile=$1}; if($2==tRNA){print $1,pdbfile}}' $PDB_dir/$pdb_lower\-chain-id-mapping.txt | awk '{print $1}'`
			echo $tRNAchainIDfinal
			tRNApdbFile=`awk -v tRNA=$tRNAchainID '{if($2==""){pdbfile=$1}; if($2==tRNA){print $1,pdbfile}}' $PDB_dir/$pdb_lower\-chain-id-mapping.txt | awk '{print $2}' | sed 's/://'`
			echo $tRNApdbFile

			#### check which ribosome the tRNA belongs to: (assuming the tRNA is in the same pdb file as its ribosome!)
			cd $Matrices_dir
			numribosomes=`ls $pdb_lower\_*.out | wc -l | awk '{print $1}'`
			if [ $numribosomes -eq 1 ]
			then
				matrixout_file=$pdb_lower\_1.out
			fi
			if [ $numribosomes -gt 1 ]
			then
				matrixout_file=`grep "$tRNApdbFile" $pdb_lower\_*.out | head -n 1 | awk '{print $1}' | sed 's/://'`
			fi
			ribosome_number=`echo "$matrixout_file" | sed "s/$pdb_lower\_//; s/.out//"`
			
			echo $matrixout_file
			echo $ribosome_number
	
			let tRNA_number=$tRNA_number+1
			NAME=$pdb_lower\_$ribosome_number\_$tRNA_number
			echo $NAME
			
			cd $SOURCE/tRNA-tcl/OUTPUT	
			sbatch $SOURCE/tRNA-tcl/drive_tRNA.bash $PDB_dir/$tRNApdbFile $tRNAchainIDfinal $NAME $matrixout_file
		done
	done
done

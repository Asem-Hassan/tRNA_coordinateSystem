#!/usr/local/bin/vmd -e 

proc Get_residue_with_largest_separation_from_CCAtip {molID chainID} {
	set tRNA_P [atomselect $molID "chain $chainID and nucleic and name P"]
	set set_of_resid [$tRNA_P get resid]
	set length_of_set_of_resid [llength $set_of_resid]
	set res_final [lindex $set_of_resid end]
	set P_final [atomselect $molID "chain $chainID and nucleic and name P and resid $res_final"]
	set P_final_xyz [join [$P_final get {x y z}]]

	set max_dist 0
	set P1_found 0

	for {set i [lindex $set_of_resid 0]} {$i<$length_of_set_of_resid} {incr i} {
		set res_i [lindex $set_of_resid $i]
		set P1 [atomselect $molID "chain $chainID and nucleic and name P and resid $res_i"]
		set P1_xyz [join [$P1 get {x y z}]]
		set diff [vecsub $P1_xyz $P_final_xyz]
		set length [veclength $diff]
		#puts "$res_i ... $length"
		if {$length > $max_dist} {set max_dist $length; set P1_found $res_i; set index_in_set_of_res $i};
	}		
	
	return [list $P1_found $index_in_set_of_res $max_dist]
}

proc length_of_tRNA {molID chainID} {
	set tRNA_P [atomselect $molID "chain $chainID and nucleic and name P"]
        set set_of_resid [$tRNA_P get resid]
	set length_of_tRNA [lindex $set_of_resid end]
        #set length_of_set_of_resid [llength $set_of_resid]   #the length function would be appropriate if the tRNA doesn't have missing residues
	return $length_of_tRNA
}

proc get_minor_axis {selection} {
        set PAs_and_evals [::Orient::calc_principalaxes $selection domass];
        set eval0 [lindex $PAs_and_evals 3 3];
        set eval1 [lindex $PAs_and_evals 3 4];
        set eval2 [lindex $PAs_and_evals 3 5];

        if {$eval0 < $eval1} {set min_eval $eval0; set min_eval_no 0};
        if {$eval1 < $eval0} {set min_eval $eval1; set min_eval_no 1};
        if {$eval2 < $min_eval} {set min_eval $eval2; set min_eval_no 2};
        switch $min_eval_no {
                0 {set minor_axis [lindex $PAs_and_evals 0]}
                1 {set minor_axis [lindex $PAs_and_evals 1]}
                2 {set minor_axis [lindex $PAs_and_evals 2]}
	};
	return $minor_axis
}

proc Assess_region_ASL {molID start end chainID} {
	#calculates the dotproduct between a given vector 
	set state "found";	
	set max_i [expr floor(($end-$start)/3)];
	
	if {$max_i > 1} {
		for {set i 0} {$i<$max_i} {incr i} {
			set st_sub [expr $start+$i]; 
			set en_sub [expr $end-$i]; 
			set subregion [atomselect $molID "chain $chainID and nucleic and resid $start to $end"];
			set minor_axis [get_minor_axis $subregion] 
			
			if {$i == 0} {set first_minor_axis $minor_axis};
			set dotpro [vecdot $minor_axis $first_minor_axis]; 
			#set dotpro [vecdot $minor_axis $minor_axis_compare]; 
			#puts $dotpro; 
			#puts "################################"
			if {abs($dotpro) < 0.985} {set state "NOTfound"; break};
		}
	} else {
		set state "NOTfound";
	}
	return $state
}

proc Assess_region_acceptor_arm {molID start1 end1 start2 end2 chainID} {
	#calculates the dotproduct between a given vector 
	set state "found";	
	set max_i [expr floor(($end2-$start2)/2)];
	
	if {$max_i > 1} {
		for {set i 0} {$i<$max_i} {incr i} {
			set st_sub1 [expr $start1]; 
			set en_sub1 [expr $end1]; 
			set st_sub2 [expr $start2+$i];
			set en_sub2 [expr $end2];
			set subregion [atomselect $molID "chain $chainID and nucleic and ((resid $st_sub1 to $en_sub1) or (resid $st_sub2 to $en_sub2))"]; 
			set minor_axis [get_minor_axis $subregion]

			if {$i == 0} {set first_minor_axis $minor_axis};
			set dotpro [vecdot $minor_axis $first_minor_axis]; 
			#puts $dotpro; 
			#puts "################################"
			if {abs($dotpro) < 0.985} {set state "NOTfound"; break};
#		}
#	} else {
#		set state "NOTfound";
#	}
	return $state
}

proc find_ASL {molID chainID} {
	set lower_res_resid [lindex [Get_residue_with_largest_separation_from_CCAtip $molID $chainID] 0]
	set lower_res_index_in_chain [lindex [Get_residue_with_largest_separation_from_CCAtip $molID $chainID] 1]
	set difference 0
	set start_found 0
	set end_found 0

	set length [length_of_tRNA $molID $chainID]	
	set tolast [expr $length - $lower_res_index_in_chain]
	set tobegin [expr $lower_res_index_in_chain]
	set maxsearch [expr min($tolast,$tobegin)]

	for {set delta 1} {$delta < $maxsearch} {incr delta} {
		set start [expr $lower_res_resid - $delta]; 
		set end [expr $lower_res_resid + $delta] 
		set state [Assess_region_ASL $molID $start $end $chainID]; 
		if {$state == "found"} {
			set diff [expr $end-$start]
			if {$diff > $difference} {
				set start_found $start
				set end_found $end
				#puts "$start ... $end ... diff=$diff"
				set difference $diff
			}
		}
	}
	return [list $start_found $end_found]

}

proc find_acceptor_arm {molID chainID end_of_ASL_proc} {
	set length [length_of_tRNA $molID $chainID]
	set max_i [expr floor($length/2)]
	set end2 [expr $length-3]
	set start_compare $end_of_ASL_proc

	for {set i 0} {$i<$max_i} {incr i} {
		set start2 [expr $end_of_ASL_proc-$i]; 
		set state [Assess_region_acceptor_arm $molID 1 6 $start2 $end2 $chainID]; 
		#puts "$start2 ... $state";
		if {$state == "found"} {
			if {$start2 < $start_compare} {
				set start_compare $start2	
			}
		}
	}
	set start_found [expr $start_compare +2];   #takes away 2 nucleotides to account for the 10 degrees threshold of error. this is now entering in the ASL.
	return [list 1 6 $start_found $end2]
}

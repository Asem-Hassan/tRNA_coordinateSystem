#!/usr/local/bin/vmd -e 

#reference tRNA coordinate system in 4v9d_2 reference structure pdb coordinates:
set X_ref {-0.843167581057584 0.37273546301225546 -0.3874747796522271}
set Y_ref {-0.4849443200701885 -0.8384277950302104 0.24873246461296983}
set Z_ref {-0.2321582147699386 0.3976268440810043 0.8876910815035914}

global X_ref Y_ref Z_ref

set pi 3.14159265358979323846 
global pi

proc invert_4by4_matrix {M} {  
	#the given input matrix M has to be a transformation matrix (4th row is 0 0 0 1)
	set M_tr [transtranspose $M]
	set T1 [lindex $M 0 3]
	set T2 [lindex $M 1 3]
	set T3 [lindex $M 2 3]
	set T_mag [veclength "$T1 $T2 $T3"]
	set delta [expr 1+$T_mag**2]
	set aux_matrix [list "1 0 0 [vecscale $T1 -1]" "0 1 0 [vecscale $T2 -1]" "0 0 1 [vecscale $T3 -1]" "[vecscale $T1 -1] [vecscale $T2 -1] [vecscale $T3 -1] $delta"]
	set M_inv [transmult $M_tr $aux_matrix]

	return $M_inv

}


proc form_4by4matrix_from_orthonormal_system {Vx Vy Vz} {
	set V1 [lappend Vx 0]
	set V2 [lappend Vy 0]
	set V3 [lappend Vz 0]
	set Va {0 0 0 1}

	set A_tr [list "$V1" "$V2" "$V3" "$Va"]; #the system's matrix transposed
	return [transtranspose $A_tr]
}

proc similarity_transform {M} {
	global X_ref Y_ref Z_ref
	set ref_system [form_4by4matrix_from_orthonormal_system $X_ref $Y_ref $Z_ref]
	set ref_system_inv [transtranspose $ref_system]

	set M_similar [transmult $ref_system_inv [transmult $M $ref_system]]
	return $M_similar
}


##################### have to transform the tRNA_coordinate_system matrix (by transforming ASL_minor_axis and acceptor_arm_minor_axis) to the 4v9d_2 reference structure pdb coordinates. ... have to use Mlarge and Mbody transformation matrices.

proc tRNA_coordinate_system {ASL_minor_axis acceptor_arm_minor_axis Mlarge Mbody} {
	global X_ref Y_ref Z_ref
	set moving_matrix [transmult [invert_4by4_matrix $Mbody] $Mlarge];  #needed to transform the vectors calculated in the model frame to the reference frame (4v9d_2) ... Mlarge and Mbody will be got from the subunit rotation script. and ASL_minor_axis and acceptor_arm_minor_axis are calculated before running the rotation script on the model (as the script moves the model). here we need Rbody inverse only which can be got from Mbody transposed.
	set ASL_minor_axis [vectrans $moving_matrix $ASL_minor_axis]
	set acceptor_arm_minor_axis [vectrans $moving_matrix $acceptor_arm_minor_axis]	
	
	#to conform with the convention of the direction of ASL_minor_axis and acceptor_arm_minor_axis:
	set ASL_dot_Z_ref [vecdot "$Z_ref" "$ASL_minor_axis"]
	set acc_dot_X_ref [vecdot "$X_ref" "$acceptor_arm_minor_axis"]
	if {$ASL_dot_Z_ref < 0} {set ASL_minor_axis [vecscale $ASL_minor_axis -1]}
	if {$acc_dot_X_ref > 0} {set acceptor_arm_minor_axis [vecscale $acceptor_arm_minor_axis -1]}
	
	set Z "$ASL_minor_axis"
	set Y [veccross "$Z" [vecscale "$acceptor_arm_minor_axis" -1]] 
	set X [veccross "$Y" "$Z"]
	
	return [form_4by4matrix_from_orthonormal_system $X $Y $Z]
}


proc calc_tRNA_orientation {tRNA_system} {
	global X_ref Y_ref Z_ref	
	global pi
	set ref_system [form_4by4matrix_from_orthonormal_system $X_ref $Y_ref $Z_ref]
	set ref_system_inv [invert_4by4_matrix $ref_system]
	set transformation [transmult $tRNA_system $ref_system_inv]
	set transformation_after_sim_transform [similarity_transform $transformation]
		
	#calculate Euler angles now from $transformation_after_sim_transform .. these calculations are based on the ACTIVE REPRESENTATION OF THE ROTATION MATRIX given in terms of Euler angles:
	set R13 [lindex $transformation_after_sim_transform 0 2]
	set R23 [lindex $transformation_after_sim_transform 1 2]
	set R33 [lindex $transformation_after_sim_transform 2 2]
	set R31 [lindex $transformation_after_sim_transform 2 0]
	set R32 [lindex $transformation_after_sim_transform 2 1]
	
	set theta [expr acos($R33)]
	if {$theta==0 || $theta==$pi} {  
		set psi 0; 
		set R11 [lindex $transformation_after_sim_transform 0 0]; 
		set R21 [lindex $transformation_after_sim_transform 1 0]; 
		set phi [expr atan2($R21,$R11)]
	} else {
		set phi [expr atan2($R13,-1*$R23)]
		set psi [expr atan2($R31,$R32)]

	#	if {$R23<0} {set phi [expr atan($R13/(-1*$R23))]}
	#	if {$R23==0} {if{$R13>0}{set phi [expr $pi/2]}; if{$R13<0}{set phi [expr 3*$pi/2]}}
	#	if {$R23>0} {set phi [expr atan($R13/(-1*$R23))+$pi]}


	#	if {$R32>0} {set psi [expr atan($R31/$R32)]}
	#	if {$R32==0} {if{$R31>0}{set psi [expr $pi/2]}; if{$R31<0}{set psi [expr 3*$pi/2]}}
	#	if {$R32<0} {set psi [expr atan($R31/$R32)+$pi]}
	}

	set phi [expr $phi*180/$pi]
	set theta [expr $theta*180/$pi]
	set psi [expr $psi*180/$pi]

	return [list $phi $theta $psi]
}

proc bring_ref_tRNA {Mlarge Mbody} {
	#this proc brings the reference tRNA from 4v9d_2 ribosome and puts it on the model for comparison with the reference tRNA. This assumes that the model has been loaded without changing its coordinates (i.e. the rotation script wasn't run on it).
	set Mlarge_inv [invert_4by4_matrix $Mlarge]
	set moving_matrix [transmult $Mlarge_inv $Mbody]
	set molID [mol new /home/asem/vmd_libraries/tRNA-tcl/ref-tRNA.pdb]
	mol modstyle 0 $molID Tube 1.300000 12.000000;   #make ref_tRNA in tube representation
	mol modcolor 0 $molID ColorID 16; 	#color ref_tRNA in black
	set ref_tRNA [atomselect $molID "all"]
	$ref_tRNA move $moving_matrix	
}

#proc apply_rotation_to_ref_tRNA {phi theta psi} {
#	#This proc takes 3 Euler angles [in degrees], forms the rotation matrix in the active representation (moving vectors), then applies this matrix to the ref_tRNA structure, for comparison the ref_tRNA should be loaded and oriented on the model by using bring_ref_tRNA procedure first.
#	global pi
#	set phi_rad [expr $phi*$pi/180]
#	set theta_rad [expr $theta*$pi/180]
#	set psi_rad [expr $psi*$pi/180]
#
#	set c_phi [expr cos($phi_rad)]
#	set s_phi [expr sin($phi_rad)]
#	set c_theta [expr cos($theta_rad)]
#	set s_theta [expr sin($theta_rad)]
#	set c_psi [expr cos($psi_rad)]
#	set s_psi [expr sin($psi_rad)]
#
#	set R00 [expr 
#}


proc main {molID chainID Mlarge Mbody} {
	lassign [find_ASL $molID $chainID] ASL_st ASL_en
	lassign [find_acceptor_arm $molID $chainID $ASL_en] acc_st1 acc_en1 acc_st2 acc_en2
	set ASL [atomselect $molID "nucleic and chain $chainID and resid $ASL_st to $ASL_en"]
	set acc [atomselect $molID "nucleic and chain $chainID and ((resid $acc_st1 to $acc_en1) or (resid $acc_st2 to $acc_en2))"]
	set ASL_minor_axis [get_minor_axis $ASL]
	set acc_minor_axis [get_minor_axis $acc]
	
	set system [tRNA_coordinate_system $ASL_minor_axis $acc_minor_axis $Mlarge $Mbody]
	lassign [calc_tRNA_orientation $system] phi theta psi
	
	#### this now represents the orientation in rotation,tilt,tilt_direction coordinates as defined in what follows. And also defines the 	      #### rotation to be in the range -180 to 180 and nearer to 0. [as phi+psi lies between -360 and 360 and also is physically congruent 	    #### to its value modulo 360.

	set tilt $theta
	set tilt_direction $phi
	set rotation [expr $phi+$psi]; #rotation is between -360 and 360
	if {$rotation < 0} {
		set rotation [expr ($rotation+360)]; #rotation is now between 0 and 360
	} 
	if {$rotation > 180} {set rotation [expr $rotation - 360]}; #rotation is now between -180 and 180.
		

	return [list $rotation $tilt $tilt_direction]
}

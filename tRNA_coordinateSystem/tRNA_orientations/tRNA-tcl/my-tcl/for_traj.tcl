#!/usr/local/bin/vmd -e 

proc calc_tRNA_orientation_traj {ASL_file acceptor_file Mlarge_file Mbody_file tRNA_name} {
	set ASL_all [open $ASL_file r]
	set acceptor_all [open $acceptor_file r]
	set Mlarge_all [open $Mlarge_file r]
	set Mbody_all [open $Mbody_file r]

	set OUTPUT [open $tRNA_name\_orientation_traj w]
	set counter 0

	while {[gets $ASL_all ASL_line] >= 0} {
		set counter [expr $counter+1]
		set ASL_minor_axis "[lindex $ASL_line 1] [lindex $ASL_line 2] [lindex $ASL_line 3]"
		gets $acceptor_all acc_line
		set acc_minor_axis "[lindex $acc_line 1] [lindex $acc_line 2] [lindex $acc_line 3]"

		gets $Mlarge_all Mlarge
		gets $Mbody_all Mbody
		
		puts $ASL_minor_axis
		puts $acc_minor_axis
		puts $Mlarge
		puts $Mbody

		
		set system [tRNA_coordinate_system $ASL_minor_axis $acc_minor_axis $Mlarge $Mbody]
		puts "SYSTEM: $system\n############################"
        	lassign [calc_tRNA_orientation $system] phi theta psi
		puts "ANGLES: $phi $theta $psi\n##################"
		
		########################################## 
		###### This redefines the tRNA orientation in terms of rotation,tilt,tilt_direction and ensures continuity of rotation: 
		set tilt $theta
		set tilt_direction $phi

		set rotation [expr $phi+$psi]; #rotation is between -360 and 360
	        if {$rotation < 0} {
	                set rotation [expr ($rotation+360)]; #rotation is now between 0 and 360
        	}
		if {$rotation > 180} {set rotation [expr $rotation - 360]}; #rotation is now between -180 and 180. this means that the only discontinuity that will be exhibited in the trajectories is when the tRNA rotation angle is around 180 and flips to -180 (highly impropable for tRNA inside the ribosome).

#################### Another method to directly check continuity with the previous rotation value:	
#		if {$counter==1} {
#			set continous_rotation $rotation
#		}
#	
#		if {$counter>1} {
#			set difference [expr $rotation - $old_rotation]
#			if {abs($difference)>270} {
#				set continous_rotation [expr $rotation+360*$old_rotation/abs($old_rotation)]
#			} else {
#				set continous_rotation $rotation
#			}
#		}
#		
#		set old_rotation $continous_rotation
#		puts $OUTPUT "$continous_rotation $tilt $tilt_direction"
#######################################################################################################	
		puts $OUTPUT "$rotation $tilt $tilt_direction"
	}
	
	close $ASL_all
	close $acceptor_all
	close $Mlarge_all
	close $Mbody_all
	close $OUTPUT


} 

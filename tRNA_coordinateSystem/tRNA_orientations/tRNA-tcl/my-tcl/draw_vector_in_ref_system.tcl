#!/usr/local/bin/vmd -e 

proc draw_vector {vector} {
	set center {-78.0 128.0 103.5}
	set scale 30
	set end [vecadd $center [vecscale 30 $vector]]
	graphics 1 line $center $end width 10
}

package require psfgen

topology top_toid_dtm1.inp
topology top_tube_former.inp

set ixmol 8
set izmol 12

set nmol [expr {$ixmol*$izmol}]

for {set sgid 1} {$sgid <= $nmol} {incr sgid} {

    set segid [concat [list "U"][list $sgid]]
    segment $segid {
        first PACE
        last CT2	
	pdb cis_planar01.pdb
    }
    coordpdb cis_planar01.pdb $segid
    guesscoord
}

writepdb tmp.pdb
writepsf tmp.psf

set zoffset 4.6
set xoffset 30.0

set xval 0.0
set yval 0.0
set zval 0.0

mol load psf tmp.psf pdb tmp.pdb 

set sgid 1


for {set kmol 1} {$kmol<=$izmol} {incr kmol} {
    for {set imol 1} {$imol<=$ixmol} {incr imol} {
	    
	    set segid [concat [list "U"][list $sgid]]
	    set mov [list $xval $yval $zval]
	    set sel [atomselect top "segid $segid"]
	    
	    if {$imol%2 == 0} {
		set com [measure center $sel weight mass] 
		set matrix [transaxis z 180]
		set matrix1 [transaxis y 180]
		$sel moveby [vecscale -1.0 $com] 
		$sel move $matrix 
		$sel move $matrix1
		$sel moveby $com 
		$sel moveby {0.0 30.0 0.0}
	   }
	   
	   $sel moveby $mov
	   incr sgid
	   set xval [expr {$xval+$xoffset}]
	    	    	    
	}
	set xval 0.0
    	set zval [expr {$zval+$zoffset}]
}

set sel [atomselect top all]
set com [measure center $sel weight mass] 
$sel moveby [vecscale -1.0 $com] 

$sel writepsf cis_planar_12_8.psf
$sel writepdb cis_planar_12_8.pdb

resetpsf

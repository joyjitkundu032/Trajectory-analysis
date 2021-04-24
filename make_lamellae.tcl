package require psfgen

topology ./../rayan_parameters/top_all22_prot-peptoid.rtf

set ixmol 3
set iymol 4
set izmol 12

set nmol [expr {$ixmol*$iymol*$izmol}]

for {set sgid 1} {$sgid <= $nmol} {incr sgid} {

    set segid [concat [list "U"][list $sgid]]
    segment $segid {
	first ACE 
	last CT2 
	pdb ./NdcNte_test_autopsf.pdb
    }
    coordpdb ./NdcNte_test_autopsf.pdb $segid
    guesscoord
}

writepdb tmp.pdb
writepsf tmp.psf

set xoffset 136.0
set yoffset 26.0
set zoffset 4.0

set xval 0.0
set yval 0.0
set zval 0.0

mol load psf tmp.psf pdb tmp.pdb 

set sgid 1


for {set kmol 1} {$kmol<=$izmol} {incr kmol} {
    for {set jmol 1} {$jmol<=$iymol} {incr jmol} {
	for {set imol 1} {$imol<=$ixmol} {incr imol} {
	    
	    set segid [concat [list "U"][list $sgid]]
	    set mov [list $xval $yval $zval]
	    set sel [atomselect top "segid $segid"]

	    if {$kmol%2 == 0} {
		set com [measure center $sel weight mass] 
		set matrix [transaxis y 180] 
		$sel moveby [vecscale -1.0 $com] 
		$sel move $matrix 
		$sel moveby $com 
		$sel moveby {59.0 14.0 0.0}
	    }

	    $sel moveby $mov
	    incr sgid
	    set xval [expr {$xval+$xoffset}]
	    	    	    
	}	

	if {$kmol%2 == 0} {
	    set xval 0.0
	} else {
	    #set xval [expr {$xoffset*0.50}]
	    set xval 0.0
	}
	set yval [expr {$yval+$yoffset}]

    }
    
    set xval 0.0
    set yval 0.0
    set zval [expr {$zval+$zoffset}]
}

set sel [atomselect top all]
set com [measure center $sel weight mass] 
$sel moveby [vecscale -1.0 $com] 

$sel writepsf lam_multi_3_4_12.psf
$sel writepdb lam_multi_3_4_12.pdb

resetpsf
#exit

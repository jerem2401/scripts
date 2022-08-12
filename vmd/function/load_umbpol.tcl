#(resname HIS or resname ARG or resname LYS) and same residue as (within 4 of nucleic)
#ions and name NA and within 4 of nucleic
#foreach mol [molinfo list] {
#mol selection nucleic
#mol addrep $mol
#mol modstyle 0 $mol NewCartoon
#mol modstyle 1 $mol NewCartoon
#mol modcolor 1 $mol ResID
#}
#vmdupol "test2_2.xtc test2.pdb ."
#vmdupol "umb2.xtc 0umb_chains.pdb E_46.840 E_47.680 E_48.520 E_49.360 E_50.200 E_50.550"
proc load_umbpol {args} {
    set mol 0
    mol new /home/jeremy/mnt/smaug/syncsim/pol/ref/oc/5iyb_clean1_correcteddna.pdb type {pdb}
    mol rename $mol oc_ref
    mol selection nucleic
    mol addrep $mol
    mol modcolor 1 $mol ResID
    mol selection (chain A and resid 320 to 338 or chain B and resid 450 to 466 or chain B and resid 476 to 500)
    mol addrep $mol
    mol modcolor 2 $mol ColorID 1
    set mol 1
    for {set x 2} {$x < [llength $args]} {incr x} {
	#loadings
        puts "loading [lindex $args $x]/[lindex $args 1] [lindex $args $x]/[lindex $args 0]"
	#mol new [lindex $args $x]/[lindex $args 1] type {pdb} first 0 last -1 step 1 waitfor all
	mol new [lindex $args $x]/[lindex $args 1] first 0 last -1 step 1 waitfor all
        mol addfile [lindex $args $x]/[lindex $args 0] type {xtc} first 1 last -1 step 1 molid $mol waitfor all
	mol rename $mol [lindex $args $x]
	#selections
        mol selection nucleic
        mol addrep $mol
        mol modcolor 1 $mol ResID

        #mol selection (protein within 30 of nucleic)
        #mol addrep $mol
	#mol selupdate 1 $mol on
        #mol modcolor 2 $mol Molecule

	mol selection (chain A and resid 320 to 338 or chain B and resid 450 to 466 or chain B and resid 476 to 500)
	mol addrep $mol
	mol modcolor 3 $mol ColorID 1
	#align first frame to ref
	set reference_sel  [atomselect 0 "((chain O and resid 76 to 86) or (chain N and resid 364 to 372) or (chain A and resid 831 to 870) or (chain H and resid 121 to 127) or (chain H and resid 140 to 146)) and backbone"]
	set comparison_sel [atomselect $mol "((chain O and resid 76 to 86) or (chain N and resid 364 to 372) or (chain A and resid 831 to 870) or (chain H and resid 121 to 127) or (chain H and resid 140 to 146)) and backbone" frame last]
	set transformation_mat [measure fit $comparison_sel $reference_sel]
	set move_sel [atomselect $mol "all" frame last]
	$move_sel move $transformation_mat
	#align all frames to last frame
	fitframes $mol "protein"
	#smoothing

	mol selection (protein and same residue as (within 4 of nucleic))
	mol representation CPK
	mol addrep $mol
	mol selupdate 3 $mol on

	mol selection (protein within 30 of nucleic)
	mol representation NewCartoon
	mol addrep $mol

	mol smoothrep $mol 0 2
	mol smoothrep $mol 1 2
	mol smoothrep $mol 2 2
	mol smoothrep $mol 3 2
	mol smoothrep $mol 4 2
	mol showrep $mol 0 off
	mol inactive $mol
	mol off $mol
	incr mol
    display resetview
    }
}
#set n [molinfo $mol get numreps]
#for {set i 0} {$i < $n} {incr i} {
#	mol showrep $mol $i off
#}

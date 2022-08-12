#foreach mol [molinfo list] {
#mol selection nucleic
#mol addrep $mol
#mol modstyle 0 $mol NewCartoon
#mol modstyle 1 $mol NewCartoon
#mol modcolor 1 $mol ResID
#}
#vmdpol "./tk7000_r1c0.4k10000_r2c1.99k20000-30000_125ns/check_stab/rep* ./tk7000_r1k10000_r2k20000-30000_125ns/check_stab/rep*"
proc load_pol {args} {
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
    foreach x $args {
	#loadings
        puts "loading $x/0.pdb $x/nopbc2.xtc"
	mol new $x/0_chains.pdb type {pdb} first 0 last -1 step 1 waitfor all
        mol addfile $x/nopbc2.xtc type {xtc} first 1 last -1 step 1 molid $mol waitfor all
	mol rename $mol $x
	#selections
        mol selection nucleic
        mol addrep $mol
        mol modcolor 1 $mol ResID

        mol selection (protein within 30 of nucleic)
        mol addrep $mol
	mol selupdate 1 $mol on
        mol modcolor 2 $mol Molecule

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

	##smoothing
	mol smoothrep $mol 0 2
	mol smoothrep $mol 1 2
	mol smoothrep $mol 2 2
	mol smoothrep $mol 3 2
	mol showrep $mol 0 off
	incr mol
    display resetview
     }
}

#foreach mol [molinfo list] {
#mol selection nucleic
#mol addrep $mol
#mol modstyle 0 $mol NewCartoon
#mol modstyle 1 $mol NewCartoon
#mol modcolor 1 $mol ResID
#}

proc load_opro {args} {
    set mol 0
    foreach x $args {
	#loadings
        puts "loading $x/0.gro $x/nopbc.xtc"
	mol new $x/0.gro type {gro} first 0 last -1 step 1 waitfor all
        mol addfile $x/nopbc.xtc type {xtc} first 1 last -1 step 1 molid $mol waitfor all
	mol rename $mol $x
	#selections
        mol selection resname FOM
        mol addrep $mol
	#mol selupdate 1 $mol on
	mol modstyle 1 $mol VDW

	fitframes $mol "protein"
	#smoothing
	mol smoothrep $mol 0 2
	mol smoothrep $mol 1 2
	incr mol
    display resetview
     }
}

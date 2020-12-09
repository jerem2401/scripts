#foreach mol [molinfo list] {
#mol selection nucleic
#mol addrep $mol
#mol modstyle 0 $mol NewCartoon
#mol modstyle 1 $mol NewCartoon
#mol modcolor 1 $mol ResID
#}

proc load_pol {args} {
    set mol 0
    foreach x $args {
        puts "loading $x/0.pdb $x/nopbc2.xtc"
	mol new $x/0_chains.pdb type {pdb} first 0 last -1 step 1 waitfor all
        mol addfile $x/nopbc2.xtc type {xtc} first 1 last -1 step 1 molid $mol waitfor all
	mol rename $mol $x
        mol selection nucleic
        mol addrep $mol
        mol modcolor 1 $mol ResID
	fitframes $mol "protein"
	mol smoothrep $mol 0 2
	mol smoothrep $mol 1 2
	incr mol
     }
}

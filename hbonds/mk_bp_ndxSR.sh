#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
#set -o errexit   # abort on nonzero exitstatus
#set -o nounset   # abort on unbound variable
#set -o pipefail  # dont hide errors within pipes

PWD=$(pwd)
if [ "$HOSTNAME" == smaug ]; then
	base=$(echo ${PWD%/simulation*})
	GMX='gmx'
	run=''
	nt='-nt 12'
else
	module load anaconda3/2020.07 && source activate env1
	base=$HOME
	source /usr/users/cmb/shared/spack/shared.bash
	module load gromacs@2021.1-mpi
	GMX='gmx_mpi'
	run='mpirun -np 12'
	nt=''
fi

traj=''
dir=''
mknd=1
ana=''

while [ $# -gt 0 ]; do
	case "$1" in
		-d)
		  shift
		  dirumb=$1
		  ;;
		-ana) shift
		  ana=$1
		  mknd=0;;
		-rtx)
		  module unload gromacs@2021.1-mpi
		  module load gromacs@2021.1
		  GMX='gmx'
		  run=''
		  nt='-nt 12'
		  ;;
	esac
	shift
done

dir="${dirumb}/short_range_2"


if (($mknd == 1)); then

	[[ -d $dir ]] && echo "$dir already exists, exiting" && exit || mkdir $dir

	for i in ${dir}; do
		cat <<- EOF > ${i}/md.mdp
		integrator              = md        ; leap-frog stochastic integrator
		;nsteps                  = 10000000
		;nsteps                  = 12500000
		nsteps                  = 43750000
		dt                      = 0.004     ; 4 fs
		; Output control
		nstxout                 = 0         ; suppress bulky .trr file by specifying 
		nstvout                 = 0         ; 0 for output frequency of nstxout,
		nstfout                 = 0         ; nstvout, and nstfout
		nstenergy               = 5000
		energygrps              = Protein DNA Water
		nstlog                  = 500000      ; update log file every 10.0 ps
		;nstxout-compressed      = 120000    ; save compressed coordinates every 10 ps
		nstxout-compressed      = 12500     ; save compressed coordinates every 10 ps
		compressed-x-grps       = ZN_MG_DNA_Protein    ; save the whole system
		; Bond parameters
		continuation            = no        ; Restarting after NPT 
		constraint_algorithm    = lincs     ; holonomic constraints 
		constraints             = h-bonds   ; bonds involving H are constrained
		lincs_iter              = 1         ; accuracy of LINCS
		lincs_order             = 4         ; also related to accuracy
		; Neighborsearching
		cutoff-scheme           = Verlet    ; Buffered neighbor searching
		ns_type                 = grid      ; search neighboring grid cells
		nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
		rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
		rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
		; Electrostatics
		coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
		pme_order               = 4         ; cubic interpolation
		fourierspacing          = 0.16      ; grid spacing for FFT
		; Temperature coupling is on
		tcoupl                  = V-rescale
		tc-grps                 = ZN_MG_DNA_Protein SOL_CL_NA       ;two coupling groups - more accurate
		tau_t                   = 0.1     0.1                       ; time constant, in ps
		ref_t                   = 300     300                       ; reference temperature, one for each group, in K
		; Pressure coupling is on
		pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
		pcoupltype              = isotropic             ; uniform scaling of box vectors
		tau_p                   = 2.0                   ; time constant, in ps
		ref_p                   = 1.0                   ; reference pressure, in bar
		compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
		; Periodic boundary conditions
		pbc                     = xyz       ; 3-D PBC
		; Dispersion correction
		DispCorr                = EnerPres  ; account for cut-off vdW scheme
		; Velocity generation
		gen_vel                 = yes        ; Velocity generation is off
		gen-seed                = -1
		EOF

		conf=$(command ls ${dirumb}/conf_*.gro)

		$GMX grompp -nice 0 -f ${i}/md.mdp -c $conf -t ${dirumb}/state.cpt -p \
		$base/simulation/syncsim/pol/heavy_h/ref/polcc_ZN_heavyH.top -n \
		$base/simulation/syncsim/pol/heavy_h/ref/index.ndx -o ${i}/srint.tpr -maxwarn 1
		wait

		$run $GMX mdrun -nice 19 -rerun ${dirumb}/traj_comp.xtc $nt -g ${i}/srint.log \
		-e ${i}/srint.edr -x ${i}/srint.xtc -s ${i}/srint.tpr
		wait

		echo 21 22 | $GMX energy -nice 0 -f ${i}/srint.edr -o ${i}/srint_prot-dna.xvg \
		&> ${i}/outnrg_prot-dna.txt
		echo 33 34 | $GMX energy -nice 0 -f ${i}/srint.edr -o ${i}/srint_dna-dna.xvg \
		&> ${i}/outnrg_dna-dna.txt
		echo 37 38 | $GMX energy -nice 0 -f ${i}/srint.edr -o ${i}/srint_dna-water.xvg \
		&> ${i}/outnrg_dna-water.txt
	done

elif [[ ! -z $ana ]]; then
	set -o noglob
	statssr.py -f $ana
	set +o noglob
fi
#for i in E*/hbond/analyzehb.txt; do a=$(dirname $i); E=$(echo ${a:0:-6}); hb=$(tail -n1 $i ); real=$(awk '{ total += $2 } END { print total/NR }' $E/colvar_${E:2}.txt); printf "%-8s %20s %20s \n" $E $real $hb >> hbout.txt; done

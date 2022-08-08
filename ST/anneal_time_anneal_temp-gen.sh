#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes

ftt=''
nsteps=''

while [ $# -gt 0 ]; do
    case "$1" in
	-temp) shift
	    temp=$1;;
	-ngroup) shift
	    ngroup=$1;;
	-dt) shift
	    dt=$1
	    echo "time-step chosen: ${dt}";;
	-ksi) shift
	    ksi=$1
	    echo "ksi = $ksi";;
	-ftemptimp) shift
	    ftt=$1
	    echo "adding ${ftt}ns to final temperature";;
	--help|-h)
	    echo "example: anneal_time_anneal_temp-gen.sh -temp '300 334 2' -ngroup 2 (-ftemptimp 100000)"
	    exit 192;;
	*) echo -e "\n$0: Error, unknown argument: $1"
	   exit 192;;
    esac
    shift
done

temp0=$(echo $temp | cut -f1 -d ' ')
tempf=$(echo $temp | cut -f2 -d ' ')
tempr=$(echo $temp | cut -f3 -d ' ')
anneal_temp=''
anneal_npnts=0

for j in $(seq 1 $ngroup); do
    for i in $(seq $temp0 $tempr $tempf); do
        anneal_temp="${anneal_temp} $i $i"
	anneal_npnts=$(($anneal_npnts+1))
    done
done

anneal_times=0
a=anneal_times
for i in $(seq 1 $((${anneal_npnts}-1))); do
    if [ $(($i%2)) -eq 0 ]; then
        b=$(($a+100))
	a=$b
	anneal_times="${anneal_times} ${b}"
    else
        b=$(($a+2000))
        a=$b
	anneal_times="${anneal_times} ${b}"
    fi
done

lasttime=$(echo $anneal_times | awk 'NF>1{print $NF}')
if [[ ! -z $ftt ]]; then
    lasttime=$(((${lasttime}+${ftt})))
    nsteps=$(echo "a=${lasttime}; b=${dt}; if ( a%b ) a/b+1 else a/b" | bc)
    trim=$(echo $anneal_times | sed '$s/\w*$//')
    anneal_times="${trim}${lasttime}"
else
    nsteps=$(echo "a=${lasttime}; b=${dt}; if ( a%b ) a/b+1 else a/b" | bc)
fi

anneal=''
tmp=''
tmp2=''
for i in $(seq 1 $ngroup); do
    anneal="${anneal} single"
    tmp="${tmp} $anneal_npnts"
    tmp2="${tmp2} ${anneal_times}"
done
anneal_npnts=$tmp
anneal_times=$tmp2

cat << EOF > anneal_options.mdp
;	                :-) GROMACS - gmx grompp, 2020.4-MODIFIED (-:
;	
;	Executable:   /home/uni05/cmb/shared/opt/gromacs/broadwell/2020.4-plumed/bin/gmx
;	Data prefix:  /home/uni05/cmb/shared/opt/gromacs/broadwell/2020.4-plumed
;	Working dir:  /home/uni08/jlapier/simulation/syncsim/opro/finalsetup/longer_Zsampling/pull_po3
;	Command line:
;	  gmx grompp -p opro.top -c npt.gro -n index.ndx -f mdout.mdp -o rep_d.z_4/md.tpr -maxwarn 1

; VARIOUS PREPROCESSING OPTIONS
; Preprocessor information: use cpp syntax.
; e.g.: -I/home/joe/doe -I/home/mary/roe
include                  = 
; e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
; define                   = -DPOSRES -DPOSRES_FC_BB=100.0 -DPOSRES_FC_SC=10.0 -DPOSRES_FC_LIPID=100 -DDIHRES -DDIHRES_FC=100.0 -DPOSRES_FOM=0

; RUN CONTROL PARAMETERS
integrator               = md-vv
; Start time and timestep in ps
tinit                    = 0
dt                       = $dt
nsteps                   = $nsteps
; For exact run continuation or redoing part of a run
init-step                = 0
; Part index is updated automatically on checkpointing (keeps files separate)
simulation-part          = 1
; mode for center of mass motion removal
comm_mode                = linear
; number of steps for center of mass motion removal
nstcomm                  = 100
; group(s) for center of mass motion removal
comm_grps                = SYSTEM

; LANGEVIN DYNAMICS OPTIONS
; Friction coefficient (amu/ps) and random seed
bd-fric                  = 0
ld-seed                  = -1

; ENERGY MINIMIZATION OPTIONS
; Force tolerance and initial step-size
emtol                    = 10
emstep                   = 0.01
; Max number of iterations in relax-shells
niter                    = 20
; Step size (ps^2) for minimization of flexible constraints
fcstep                   = 0
; Frequency of steepest descents steps when doing CG
nstcgsteep               = 1000
nbfgscorr                = 10

; TEST PARTICLE INSERTION OPTIONS
rtpi                     = 0.05

; OUTPUT CONTROL OPTIONS
; Output frequency for coords (x), velocities (v) and forces (f)
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
; Output frequency for energies to log file and energy file
nstlog                   = 2000000
nstcalcenergy            = 100
nstenergy                = 1000
; Output frequency and precision for .xtc file
nstxout-compressed       = 25000
compressed-x-precision   = 1000
; This selects the subset of atoms for the compressed
; trajectory file. You can select multiple groups. By
; default, all atoms will be written.
compressed-x-grps        = 
; Selection of energy groups
energygrps               = 

; NEIGHBORSEARCHING PARAMETERS
; cut-off scheme (Verlet: particle based cut-offs)
cutoff-scheme            = Verlet
; nblist update frequency
nstlist                  = 20
; Periodic boundary conditions: xyz, no, xy
pbc                      = xyz
periodic-molecules       = no
; Allowed energy error due to the Verlet buffer in kJ/mol/ps per atom,
; a value of -1 means: use rlist
verlet-buffer-tolerance  = 0.005
; nblist cut-off        
rlist                    = 1.2
; long-range cut-off for switched potentials

; OPTIONS FOR ELECTROSTATICS AND VDW
; Method for doing electrostatics
coulombtype              = pme
coulomb-modifier         = Potential-shift-Verlet
rcoulomb-switch          = 0
rcoulomb                 = 1.2
; Relative dielectric constant for the medium and the reaction field
epsilon-r                = 1
epsilon-rf               = 0
; Method for doing Van der Waals
vdwtype                  = Cut-off
vdw-modifier             = Force-switch
; cut-off lengths       
rvdw_switch              = 1.0
rvdw                     = 1.2
; Apply long range dispersion corrections for Energy and Pressure
DispCorr                 = No
; Extension of the potential lookup tables beyond the cut-off
table-extension          = 1
; Separate tables between energy group pairs
energygrp-table          = 
; Spacing for the PME/PPPM FFT grid
fourierspacing           = 0.12
; FFT grid size, when a value is 0 fourierspacing will be used
fourier-nx               = 0
fourier-ny               = 0
fourier-nz               = 0
; EWALD/PME/PPPM parameters
pme-order                = 4
ewald-rtol               = 1e-05
ewald-rtol-lj            = 0.001
lj-pme-comb-rule         = Geometric
ewald-geometry           = 3d
epsilon-surface          = 0
implicit-solvent         = no

; OPTIONS FOR WEAK COUPLING ALGORITHMS
; Temperature coupling  
tcoupl                   = v-rescale
nsttcouple               = -1
nh-chain-length          = 10
print-nose-hoover-chain-variables = no
; Groups to couple separately
tc_grps                  = Protein_POPE Water_and_ions_FOM
; Time constant (ps) and reference temperature (K)
tau_t                    = 1.6 1.6
ref_t                    = 300 300
; pressure coupling     
pcoupl                   = Parrinello-Rahman
pcoupltype               = semiisotropic
nstpcouple               = -1
; Time constant (ps), compressibility (1/bar) and reference P (bar)
tau_p                    = 5.0
compressibility          = 4.5e-5  4.5e-5
ref_p                    = 1.0     1.0
; Scaling of reference coordinates, No, All or COM
refcoord_scaling         = com

; OPTIONS FOR QMMM calculations
QMMM                     = no
; Groups treated Quantum Mechanically
QMMM-grps                = 
; QM method             
QMmethod                 = 
; QMMM scheme           
QMMMscheme               = normal
; QM basisset           
QMbasis                  = 
; QM charge             
QMcharge                 = 
; QM multiplicity       
QMmult                   = 
; Surface Hopping       
SH                       = 
; CAS space options     
CASorbitals              = 
CASelectrons             = 
SAon                     = 
SAoff                    = 
SAsteps                  = 
; Scale factor for MM charges
MMChargeScaleFactor      = 1

; GENERATE VELOCITIES FOR STARTUP RUN
gen-vel                  = yes
gen-temp                 = 300
gen-seed                 = -1

; OPTIONS FOR BONDS    
constraints              = h-bonds
; Type of constraint algorithm
constraint_algorithm     = LINCS
; Do not constrain the start configuration
; continuation             = no
; Use successive overrelaxation to reduce the number of shake iterations
Shake-SOR                = no
; Relative tolerance of shake
shake-tol                = 0.0001
; Highest order in the expansion of the constraint coupling matrix
lincs-order              = 4
; Number of iterations in the final step of LINCS. 1 is fine for
; normal simulations, but use 2 to conserve energy in NVE runs.
; For energy minimization with constraints it should be 4 to 8.
lincs-iter               = 1
; Lincs will write a warning to the stderr if in one step a bond
; rotates over more degrees than
lincs-warnangle          = 30
; Convert harmonic bonds to morse potentials
morse                    = no

; ENERGY GROUP EXCLUSIONS
; Pairs of energy groups for which all non-bonded interactions are excluded
energygrp-excl           = 

; WALLS                
; Number of walls, type, atom types, densities and box-z scale factor for Ewald
nwall                    = 0
wall-type                = 9-3
wall-r-linpot            = -1
wall-atomtype            = 
wall-density             = 
wall-ewald-zfac          = 3

; COM PULLING          
pull                     = yes
pull_nstxout             = 200
pull_nstfout             = 200
pull_ncoords             = 3
pull_ngroups             = 4

pull_group1_name         = custom_calpha
pull-group1-pbcatom      = 2146
pull_group2_name         = cognoh
pull-group2-pbcatom      = 94397 
pull-pbc-ref-prev-step-com = yes

pull_coord1_type         = flat-bottom
pull_coord1_geometry     = distance
pull_coord1_dim          = Y Y N
pull_group3_name         = cogfom
pull-group3-pbcatom      = 94397 
pull_coord1_groups       = 1 3
pull_coord1_start        = no
pull_coord1_init         = 1.0
pull_coord1_rate         = 0.0
pull_coord1_k            = 1000

pull_coord2_type        = flat-bottom
pull_coord2_geometry    = angle-axis
pull_coord2_dim         = Y Y Y
pull_group4_name        = cogpo3
pull-group4-pbcatom     = 94397 
pull_coord2_vec         = 0 0 1
pull_coord2_groups      = 2 4
pull_coord2_start       = no
pull_coord2_init        = 45
pull_coord2_rate        = 0.0
pull_coord2_k           = 7878.7

pull_coord3_type        = umbrella
pull_coord3_geometry    = direction
pull_coord3_dim         = N N Y
pull_coord3_vec         = 0.0 0.0 -1.0
pull_coord3_groups      = 1 2
pull_coord3_start       = no
pull_coord3_init        = $ksi
pull_coord3_rate        = 0.0
pull_coord3_k           = 2000


; AWH biasing          
awh                      = no

; ENFORCED ROTATION    
; Enforced rotation: No or Yes
rotation                 = no

; Group to display and/or manipulate in interactive MD session
IMD-group                = 

; NMR refinement stuff 
; Distance restraints type: No, Simple or Ensemble
disre                    = No
; Force weighting of pairs in one distance restraint: Conservative or Equal
disre-weighting          = Conservative
; Use sqrt of the time averaged times the instantaneous violation
disre-mixed              = no
disre-fc                 = 1000
disre-tau                = 0
; Output frequency for pair distances to energy file
nstdisreout              = 100
; Orientation restraints: No or Yes
orire                    = no
; Orientation restraints force constant and tau for time averaging
orire-fc                 = 0
orire-tau                = 0
orire-fitgrp             = 
; Output frequency for trace(SD) and S to energy file
nstorireout              = 100

; Free energy variables
free-energy              = no
couple-moltype           = 
couple-lambda0           = vdw-q
couple-lambda1           = vdw-q
couple-intramol          = no
init-lambda              = -1
init-lambda-state        = -1
delta-lambda             = 0
nstdhdl                  = 50
fep-lambdas              = 
mass-lambdas             = 
coul-lambdas             = 
vdw-lambdas              = 
bonded-lambdas           = 
restraint-lambdas        = 
temperature-lambdas      = 
calc-lambda-neighbors    = 1
init-lambda-weights      = 
dhdl-print-energy        = no
sc-alpha                 = 0
sc-power                 = 1
sc-r-power               = 6
sc-sigma                 = 0.3
sc-coul                  = no
separate-dhdl-file       = yes
dhdl-derivatives         = yes
dh_hist_size             = 0
dh_hist_spacing          = 0.1

; Non-equilibrium MD stuff
acc-grps                 = 
accelerate               = 
freezegrps               = 
freezedim                = 
cos-acceleration         = 0
deform                   = 

; simulated tempering variables
simulated-tempering      = no
simulated-tempering-scaling = geometric
sim-temp-low             = 300
sim-temp-high            = 300

; Ion/water position swapping for computational electrophysiology setups
; Swap positions along direction: no, X, Y, Z
swapcoords               = no
adress                   = no

; User defined thingies
user1-grps               = 
user2-grps               = 
userint1                 = 0
userint2                 = 0
userint3                 = 0
userint4                 = 0
userreal1                = 0
userreal2                = 0
userreal3                = 0
userreal4                = 0
; Electric fields
; Format for electric-field-x, etc. is: four real variables:
; amplitude (V/nm), frequency omega (1/ps), time for the pulse peak (ps),
; and sigma (ps) width of the pulse. Omega = 0 means static field,
; sigma = 0 means no pulse, leaving the field to be a cosine function.
electric-field-x         = 0 0 0 0
electric-field-y         = 0 0 0 0
electric-field-z         = 0 0 0 0

; Density guided simulation
density-guided-simulation-active = false

; Simulated Annealing
annealing                                = ${anneal}
annealing_npoints                        = ${anneal_npnts}
annealing_time                           = ${anneal_times}
annealing_temp                           = ${anneal_temp}
EOF


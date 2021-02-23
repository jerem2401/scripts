#!/bin/bash
shopt -s -o nounset
function =() { awk 'BEGIN {print '"$*"'}' < /dev/null; }

do_submit()
{
    local jobfile=$1
    case $Qsystem in
	LSF)
	    bsub < $jobfile
	    ;;
	moab|PBS)
	    msub $jobfile
	    ;;
	slurm)
	    sbatch $jobfile
	    ;;
	*)
	    echo "Unkonwn Qsystem = $Qsystem"; exit 1
	    ;;
    esac
}


STARTDIR=`pwd`
MDRUN=mdrun
GO=""
DEFFNM=""
OPT=""
jname=md
walltime='24:00:00'
w=""
nnodes=1
rc=""
tune=""
ppn=""
ntArg=""
verb=""
extraLine=""
nperjob=0
nchain=0
GOchain=""
extend=0
opt=""

while [ $# -gt 0 ]; do
    case "$1" in
	-pat) shift
	    dir="$1" ;;
	-nperjob) shift
	    nperjob=$1 ;;
	-nchain) shift
	    nchain=$1 ;;
	-ext)
	    extend=1;;
	*)
	    # options may contain multiple words -
	    # pack into one word as a work-around to avoid trouble in jobscript.sh	   
	    opt="$opt $(echo $1| sed 's/ /@/g')"
	esac
	shift
done


# find out queueing systems
Qsystem=unset
if hostname | grep -q gwdu; then
    echo "This is at the GWDG"
    Qsystem=slurm
    if [ $nchain -gt 1 ]; then
        jobscript="/usr/users/jlapier/gitrepo/scripts/jobscript/jobchain_slurm_gwdg.sh -n $nchain"
        bChain=1
    else
        bChain=0
        jobscript=/usr/users/jlapier/gitrepo/scripts/jobscript/jobscript_slurm_gwdg.sh
    fi
elif hostname | egrep -q 'jj28l..'; then
    echo This is JUROPA
    Qsystem=moab
elif hostname | egrep -q 'hicegate'; then
    Qsystem=PBS
elif [ `hostname` = "smaug" ]; then
    echo "This is on smaug"
    Qsystem=slurm
    if [ $nchain -gt 1 ]; then
        jobscript="/home/users/jeremy/gitrepo/scripts/jobscript/jobchain_slurm.sh -n $nchain"
        bChain=1
    else
        bChain=0
        #jobscript=/data/users/jeremy/gitrepo/scripts/jobscript/common_jobscript_smaug.sh
        jobscript=/home/users/jeremy/gitrepo/scripts/jobscript/jobscript_slurm_smaug.sh
    fi
else
    echo "$0: ERROR. Don't know this machine ($(hostname))"; exit 1
fi
echo "Found queueing system $Qsystem"

N_RUNS=0
dirlist=""
for i in `ls -d $dir 2> /dev/null`; do
    # Check if run is alredy finished:
    if [  -e $i/confout.gro  ] && (($extend == 0)); then 
	echo "Found $i/confout.gro. Not running this directory."
	continue
    else
	let N_RUNS++
	dirlist="$dirlist $i"
    fi
done

N_RUNS=`ls -d $dir 2> /dev/null | wc -l`
echo "Preparing jobscripts for $N_RUNS runs"

if [ "$nperjob" != 0 ]; then
    njobs=$[(N_RUNS+nperjob-1)/nperjob]
    echo "Running $nperjob runs within each job. This will require $njobs jobs"
    rest=$[N_RUNS%nperjob]
    if [ $rest != 0 ]; then
	nfree=$[nperjob-rest]
	echo -e "\n**********************************************************\nWARNING !!!\n"
	echo -e "CPU cores for $nfree runs are not used. You are wasting CPU recources !!"
	echo -e "*********************************************************\n"
    fi
fi

ls -d $dirlist

if [ "$nperjob" = 0 ]; then
    # We have one job (and jobscript) per directory
    for i in $dirlist; do
	cd $i
	# extend jobname (behind -j) by .$i:
	thisopt=$(echo $opt | awk '{for (i=1;i<=NF; i++) if ($i == "-j"){$(i+1)=$(i+1)".'$i'"}; print}')

	$jobscript $thisopt || exit 1
	cd $STARTDIR
    done
else
    # We run multiple runs within one job - note that this is in principle not required at the GWDG. Unfortunately, the LSF
    # does not pin threads, so we need to do that
    # turn into array
    lista=( $dirlist )
    echo
    for ((j=0; j<njobs; j++)); do
	i0=$[j*nperjob]
	nthis=$nperjob
	echo nthis = $nthis
	if [ $[i0+nthis] -gt $N_RUNS ]; then
	    nthis=$[N_RUNS-i0]
	fi
	echo Now: nthis = $nthis
	thislist=$(echo "${lista[@]:i0:nthis}" | sed 's/ /@/g')
	echo thislist = $thislist
	# extend jobname (behind -j) by .$j:       
	thisopt=$(echo $opt | awk '{for (i=1;i<=NF; i++) if ($i == "-j"){$(i+1)=$(i+1)".'$j'"}; print}')
	[ $bChain = 1 ] && ext="-ext $j" || ext=""

	# jobscript.sh OR jobchian.sh -n nchain
	echo "Writing jobfile nr $j for these dirs: $thislist"
	$jobscript $ext $thisopt -bEnsureFullNode -multidir "$thislist" || exit 1
    done
fi

exit 0




if [[ ( "$GOchain" = -go ) && ( $nchain -gt 1 ) ]]; then
    echo -e "\nJobs were submitted within jobchain.sh\n"
    exit 0
fi

if [ "$GO" != "go"  ]; then 
    echo -e "\nExiting after creating run scripts..."
    exit 0
fi

# submit everything
cd $STARTDIR
{
    if [ $nperjob = 0 ]; then
	for i in $dirlist; do
	    cd $i
	    echo "dir $i, submitting jobscript job.sh"
	    do_submit job.sh
	    cd $STARTDIR
	done
    else
	for ((j=0; j<njobs; j++)); do
	    jobfile=job.$jname.$j.sh
	    do_submit $jobfile
	done
    fi
} | tee $STARTDIR/job-ids.log



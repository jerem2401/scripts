#!/bin/bash

source ~all-jh/bin/functions.bash

shopt -s -o nounset


calc() { awk 'BEGIN {print '"$*"'}' < /dev/null; }

backoff() {
local n=0
if [ -f $1 ]
then
    while [ -f \#$1.$n\# ]
    do
        ((n++))
    done
    mv $1 \#$1.$n\#
    echo "Back off! Renamed $1 to #$1.$n#!"
fi
}

update=''
version=''
jobname=cmbjob
tpr=""
npme=""
dd=""
deffnm=""
nnodes=1
ppn=6    # Default for fang
md=1
exe=""
opt=""
mpi=0
gmxrc=gmxrc_not_set
#mdrun=/usr/users/cmb/shared/bin/run_mdrun.sh
mdrun='gmx mdrun'
bGaussian=0
verb="-v"
bGo=0
key=""
email=""
line=""
modules=""
initmpi=""
queue="not-set"
launch=""
bGMXRC=0
machine=XX
bMPI=0
qExtension=""
depline=""
gpu_shares=""
gpu_id=""
Qsystem="none"
stepout=5000
feature=ice1
dir=`pwd`
multidirs=""
bPin=0
pinstride=1
nt=""
nameext=""
ptile=unset
np=""
ncores=""
bEmpty=0
bLatest=0
excludeNodes=''
bEnsureFullNode=0  # if you run 8 2-core jobs (with pinning, make sure the node is filled by your jobs): #BSUB -R np16
batchInitLine=''
gmxrcLine=''
logicalCoresPerPhysical=2
niceLevel=0
pbcfix=0

sbatch_tempfile=`mktemp sbatch.tempXXXXX`
#rm -f $sbatch_tempfile
trap "rm -rf $sbatch_tempfile" EXIT


case `hostname` in
    smaug|fang?|fang??)
        ppn=6
        gmxrc=/data/shared/opt/gromacs/2020.4/bin/GMXRC.bash
        queue=deflt
        machine=smaug
        ;;
    *)
        echo -e "WARNING: Default ppn and gmxrc unknown (hostname = `hostname`).\n"
        exit 1
esac


#Set defaults
case $machine in
    smaug|sandybridge)
        Qsystem=slurm
        walltime='2-00:00'
        maxh='-maxh 48'
        deplineTempl='#SBATCH --depend=afterok:JOBID'
        Qkey='#SBATCH'
        ;;
    *)
        echo "Unknown machine $machine"; exit 1
esac

echo "Found machine = $machine -- Qsystem $Qsystem, default walltime $walltime"


mdrun_line=""
ntFlag=""
bNoNtFlag=0
bCpi=1
scratch=0
bReqScratch=0
gpuGeneration='pascal'
plumed=""
ngpu='1'

while [ $# -gt 0 ]; do
    case "$1" in
    -empty) bEmpty=1;;
    -go) bGo=1;;
    -p) shift
            queue=$1 ;;
    -nocpi)
        bCpi=0 ;;
    -noverb)
        verb="" ;;
    -mpi)
        bMPI=1 ;;
    -no-nt-flag)
        bNoNtFlag=1 ;;
    -exe) shift
        exe=$1
        md=0
        ;;
    --time|-t)
        shift
        d=$(calc "int(($1)/24)")
        h=$(calc "int(($1)%24)")
        m1=$(calc "(($1)%24)-$h")
        m2=$(calc "int($m1*60)")
        if [ $Qsystem = slurm ]; then
        #slurm want D-HH:MM
        if [ $h -lt 10 ] & [ $m2 -lt 10 ]; then
            walltime="$d-0$h:0$m2"
        elif [ $h -lt 10 ] & [ $m2 -ge 10 ]; then
            walltime="$d-0$h:$m2"
        elif [$h -ge 10] & [ $m2 -lt 10 ]; then
            walltime="$d-$h:0$m2"
        else
            walltime="$d-$h:$m2"
        fi
        else
        echo "Unknown Qsystem = $Qsystem"; exit 1
        fi
        maxh="-maxh $1"
        echo walltime=$walltime
        echo maxh="$maxh"
        ;;
    -nice)
	shift
	niceLevel=$1 ;;
    -jobname|-J|-j)
        shift
        jobname=$1
        echo jobname=$1
        ;;
    -tpr)
        shift
        tpr=$1
        echo $tpr=$1
        ;;
    -N|--nodes|-nnodes)
        shift
        nnodes=$1
        echo nnodes=$nnodes
              #if [[ $2 =~ ^[0:9]* ]]; then
        #      shift
        #   maxnodes=$1
        #   echo maxnodes=$1
        #else
        #    echo "No range of nodes given. Using $nnodes."
        #fi
        ;;
    -ppn)
        shift
        ppn=$1
        echo ppn=$ppn
        ;;
    -bEnsureFullNode|-full)
            bEnsureFullNode=1;;
    -npme )
            shift
            npme="-npme $1"
            echo npme=$npme
            ;;
        -dd )
            shift
            dd="-dd $1"
            echo dd=$dd
            ;;
        -deffnm )
            shift
            deffnm="-deffnm $1"
            echo deffnm=$1
            ;;
        -m ) shift
             # unpack options into several words, if packed into one word with '@'
             opt=$(echo "$1" | sed 's/@/ /g')
             echo "Found extra mdrun options: opt = $opt"
             ;;
        -rc) shift
             gmxrc=$1
             ;;
        -mdrun) shift
                mdrun=$(echo "$1" | sed 's/@/ /g')
                ;;
        -mdrun_line)
	    shift
            mdrun_line=$(echo "$1" | sed 's/@/ /g')
            ;;
        -nt) shift
             nt=$1
             ntFlag="-nt $1" ;;
        # Note: this is the number of MPI processes. May be smaller than ppn when using
        # multiple OpenMP theads per MPI process
        -np) shift
             np=$1;;
        # May have to be set when using OpenMP
        -ncores)
            shift
            ncores=$1 ;;
        -opt )
            shift
            bOptimize=1
            ;;
        -stepout)
            shift
            stepout=$1 ;;
        -key)
            shift
            key=$1 ;;
        -line)
            shift
            line=$(echo "$1" | sed 's/@/ /g') ;;
        -latest)
            bLatest=1 ;;
        -gpu-generation|-gpu-gener|-g)
            shift
            gpuGeneration=$1 ;;
        -gpu_id) shift
                 gpu_id=$1 ;;
        -qshort)
            qExtension=" --qos=short" ;;
        -qlong)
            qExtension=" --qos=long" ;;
        -launch)
            launch="-launch" ;;
        -dep) shift
              depline=$(echo "$deplineTempl" | sed 's/JOBID/'$1'/g')
              echo "Using depline = $depline"
              ;;
        -email)
            shift
            email="$1" ;;
        -feature)
            shift
            feature=$1 ;;
        -multidir)
            shift
            multidirs=$(echo "$1" | sed 's/@/ /g') ;;
        -ptile)
            shift
            ptile=$1 ;;
        -pin)
            bPin=1 ;;
	-pinstride)
	    shift
	    pinstride=$1 ;;
        -scratch)
            scratch=1 ;;
        -req-scratch)
            bReqScratch=1 ;;
        -exclude-nodes)
            shift
            excludeNodes="$(echo "$1" | sed 's/@/ /g')" ;;
	-plumed)
	    shift
	    plumed="-plumed $1" ;;
	-version)
	    shift
            version=$1
	    if [ "$version" = "2018.8_plumed" ]; then
		gmxrc="/data/shared/opt/gromacs/$version/gromacs-2018.8/bin/GMXRC"
	    else
		gmxrc="/data/shared/opt/gromacs/$version/bin/GMXRC"
	    fi
            ;;
        -batch)
            shift
	    exnode=$1
            batchInitLine=$(echo "#SBATCH $exnode")
            ;;
        -update_gpu)
	    update="-update gpu";;
        -ngpu)
            ngpu=$1 ;;
    	-pbcfix)
            pbcfix=1 ;;
        *)
            echo -e "\n$0: Error, unknown argument: $1"
            exit 192
            ;;
    esac
    shift
done

if [ "$np" = "" ]; then
    np=$[nnodes*ppn]
fi

echo verb="$verb"
gmxrcLine="source $gmxrc"
echo "Using GROMACS: $gmxrc"
# pick tpr file
if [ "$multidirs" = "" ]; then
    if [[ "$tpr" = "" && "$exe" = "" && $bEmpty = 0 && "$mdrun_line" = "" ]]; then
        if [ $(ls *.tpr 2> /dev/null |wc -l) != 1  ]; then
            echo No tpr or more than one found. See: >&2
            ls *.tpr >&2
            exit 1
        else
            tpr=$(ls *tpr)
        fi
    fi
fi

# pick plumed.dat file
if [ "$multidirs" = "" ] && [ "$plumed" = "-plumed on" ]; then
    if [[ "$exe" = "" && "$mdrun_line" = "" ]]; then
        if [ $(ls *.dat 2> /dev/null |wc -l) != 1  ]; then
            echo No plumed.dat or more than one found. See: >&2
            ls *.dat >&2
            exit 1
        else
            plumed="-plumed $(ls *.dat)"
        fi
    fi
fi

nedi=$(ls *edi 2>/dev/null | wc -l)
if [ $nedi -eq 1 ]; then
    edi="-ei $(ls *edi)"
elif [ $nedi -eq 0 ]; then
    edi=""
else
    echo "There are $nedi edi files in current dir. ??" >&2
fi

ldPath=''

case $machine in
    smaug)
        [ $ptile = unset ] && ptile=$ppn
        spanline="#SBATCH --tasks-per-node=[$ptile]"
        if [ $bMPI = 1 ]; then
            mpirun="mpirun  -np $np "
            ldPath="$ldPath:"
            mdrun="${mdrun}_mpi"
        else
            mpirun=""
        fi
        ;;
           *)
        echo -e "\nError: Unknown machine $machine"; exit 1
        ;;
esac

initLD='export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$ldPath"

[ $bCpi = 1 ] && cpiArg="-cpi" || cpiArg=""

if [ "$multidirs" = "" ] ;then
    jobfile=job${key}.sh
else
    jobfile=job.${jobname}${key}.sh
fi
backoff $jobfile
############################################################################################################################

if [ "$gpuGeneration" = pascal ]; then
    excludeLine="#SBATCH --exclude=fang[41-50]"
    nGPUsPerNode=1
elif [ "$gpuGeneration" = turing ]; then
    excludeLine="#SBATCH --exclude=fang[1-40],fang50"
    nGPUsPerNode=4
elif [[ ( "$gpuGeneration" = "" ) || ( "$gpuGeneration" = any ) ]]; then
    excludeLine=""
else
    echo "ERROR, invalid GPU generation: $gpuGeneration"; exit 1
fi


############################################################################################################################

if [[ "$Qsystem" = slurm ]]; then
    if [ $bEnsureFullNode = 1 ]; then
        if [[ ( ($ncores = 12) || ( $ncores = 16 ) || ( $ncores = 20 ) || ( $ncores = 24 ) ) && ( $queue = mpi ) ]]; then
            echo "#SBATCH --batch=\"ncpus=${ncores}\"" >> $sbatch_tempfile
            echo "Will make sure that the node is full, using purely nodes with ${ncores} cores"
        fi
    fi
    {
    cat <<EOF ##SBATCH -N $nnodes temporarly removed, because of strange behaviour of the queueing system
#!/bin/bash
#SBATCH -p $queue$qExtension
#SBATCH -o $dir/myjob${key}.out
#SBATCH -e $dir/myjob${key}.err
#SBATCH -c $(echo $logicalCoresPerPhysical*$ppn | bc)
#SBATCH -t $walltime
#SBATCH --job-name=$jobname$key
#SBATCH --mail-user=$email
#SBATCH -N $nnodes
#SBATCH --nice=$niceLevel
#SBATCH --gpus=$ngpu
$excludeLine

$batchInitLine
$depline

EOF

cat $sbatch_tempfile

cat <<EOF
#echo Starting job at \$(date)
echo Jobid \$SLURM_JOBID
echo Host \$SLURM_JOB_NODELIST
echo Jobname \$SLURM_JOB_NAME
echo Subcwd \$SLURM_SUBMIT_DIR

cd $dir
$gmxrcLine
$line
$initLD
EOF
    } > $jobfile

        # exclude certain nodes (e.g. if they are known to cause trouble)
    #for i in $excludeNodes; do
     #   echo "#BSUB -R \"select[hname!='$i']\"" >> $jobfile
    #done

    #SLURM_JOB_NODELIST for SLURM_HOSTS from https://doc.itc.rwth-aachen.de/display/CC/Slurm+environmental+variables
else
    echo "Unknown queuing system = $Qsystem"
    exit 1
fi

if [ $bMPI = 0 ]; then
    if [[ ( "$ntFlag" = "" ) && ( "$Qsystem" = slurm ) ]]; then
        ntFlag="-nt \$[SLURM_CPUS_PER_TASK]"
        {
            echo "echo SLURM_CPUS_PER_TASK = \$SLURM_CPUS_PER_TASK"
        } >> $jobfile
    fi
fi

if [ $bEmpty = 1 ]; then
    echo "Wrote jobscript without gmx mdrun (empty)"
    exit 0
fi

{

    if [ "$mdrun_line" != "" ]; then

        echo "Found complete mdrun command line (option -mdrun_line), do not write mdrun line automatically" >&2
        echo "$mdrun_line &> md$key.lis"

    else

        pinArgs=""
        if [ "$multidirs" = "" ]; then
            [ "$bPin" = 1 ] && pinArgs="-pin on -pinoffset 0 -pinstride $pinstride"
            if [[ "$md" = 1 ]]; then
                #if [[ ( "$gpu_id" != unset ) ]]; then
                #    [ "$gpu_id" = '' ] && gpuID_flag="-gpu_id 0" || gpuID_flag="-gpu_id $gpu_id"
                #else
                #    gpuID_flag=""
                #fi
		gpuID_flag=""

                mdrunCall="$mpirun$mdrun"
		if ((pbcfix == 1)); then
		    maxhh=$(echo $maxh | grep -o "[0-9]*")
		    maxhhm=$(echo $maxhh - 0.05 | bc)
		    maxh="-maxh $maxhhm"
		    pldgiven=$(echo $plumed | grep -oP '(?<=\-plumed )plumed_.*\.dat')
		    if [ "$pldgiven" == plumed_nopbc.dat ]; then
		        dec="declare -a ref=( plumed_nopbc.dat plumed_pbc.dat plumed_whole.dat )"
		    elif [ "$pldgiven" == plumed_pbc.dat ]; then
		        dec="declare -a ref=( plumed_pbc.dat plumed_whole.dat plumed_nopbc.dat )"
		    else
		        dec="declare -a ref=( plumed_whole.dat plumed_pbc.dat plumed_nopbc.dat )"
		    fi
		    read -r -d '' pbcfixtxt <<EOM

START_TIME=\$(date +%s)
t2date=\$(date -d @\$START_TIME)
echo "\$t2date" >> pbcfix.log

if [ -f ./md.log ]; then
	lstec=\$(grep -oP '(?<=gmx was: ).*' pbcfix.log | tail -1)
	if (("\$lstec" == 0)); then
		pld=\$(grep -o '\-plumed plumed_.*\.dat' md.log | tail -1)
	else
		pld='$plumed'
	fi
else
	pld='$plumed'
fi

echo "1st pld used in chain is: \$pld" >> pbcfix.log

$mdrunCall $ntFlag $cpiArg -stepout $stepout $verb -s $tpr $maxh $dd $npme $deffnm $opt $edi $pinArgs \$pld $update $gpuID_flag >& md$key.lis

ecode=\$(echo \$?)
echo "exit code of 1st gmx was: \$ecode" >> pbcfix.log
$dec
iter=1

while (("\$ecode" != 0)) && (("\$iter" <= 5)); do
	((iter+=1))
	readarray -t usedpld <<< \$(grep -oP '(?<=\-plumed )plumed_.*\.dat' md.log)
	lastpld="\${usedpld[-1]}"
	echo "lastpld was \$lastpld" >> pbcfix.log
	for idx in "\${!ref[@]}"; do
		if [ "\$lastpld" == "\${ref[\$idx]}" ]; then
			nidx=\$((idx+1))
		fi
	done
	if ((\$nidx == 3)); then
		nidx=0
	fi
	nextpld="\${ref[\$nidx]}"
	echo "nextpld is \$nextpld" >> pbcfix.log

	END_TIME=\$(date +%s)
	t2date=\$(date -d @\$END_TIME)
	echo "\$t2date" >> pbcfix.log
	ELAPSED=\$(echo "scale=2; (\$END_TIME - \$START_TIME)/60/60" | bc)
	if [ \$(echo "\$ELAPSED<0.05"| bc) -eq 1 ]; then
		ELAPSED=0.05
	fi
	maxhhm=\$(echo "$maxhhm - \$ELAPSED" | bc)
	maxh="-maxh \$maxhhm"
	wait
	$mdrunCall $ntFlag $cpiArg -stepout $stepout $verb -s $tpr \$maxh $dd $npme $deffnm $opt $edi $pinArgs -plumed \$nextpld $update $gpuID_flag >& md$key.lis
	ecode=\$(echo \$?)
	echo "exit code of \${iter}th gmx was: \$ecode" >> pbcfix.log
done
EOM
		echo "$pbcfixtxt"
		else
                    mdrunArgs="$cpiArg -stepout $stepout $verb -s $tpr $maxh $dd $npme $deffnm $opt $edi $pinArgs $plumed $update"
                    echo "$mdrunCall $ntFlag $mdrunArgs $gpuID_flag >& md$key.lis"
		fi
            else
                echo "$mpirun" "$exe"
            fi
        else
	    #
	    # Running multiple dirs in one job
	    #
            idir=0
	    ndir=$(echo $multidirs | wc -w)
	    echo "Found multidirs = $multidirs" >&2
            for sdir in $multidirs; do
                if [ ! -d $dir/$sdir ]; then
                    echo "No such directory: $dir/$sdir" >&2; exit 1
                fi
                cd $dir/$sdir
                echo "cd $dir/$sdir"
                ntpr=$(ls *.tpr 2>/dev/null | wc -l)
                if [ $ntpr -ne 1 ]; then
                    echo "Found $ntpr tpr files in dir `pwd`" >&2; exit 1
                fi
                tpr=$(ls *.tpr)
                if [ "$bPin" = 1 ]; then
                    pinArgs="-pin on -pinoffset $[idir*nt*pinstride] -pinstride $pinstride"
                fi

                # Get number of mdruns running per GPU. Round up, important if ndir is an odd number
                nDirPerGPU=$(= "$ndir/$nGPUsPerNode+0.01" | round_ndig 0)
                # nDirPerGPU=$[ndir/nGPUsPerNode]
                gpuID_flag="-gpu_id $[idir/nDirPerGPU]"
                
		if [ "$plumed" != "" ]; then
		    npld=$(ls *.dat 2>/dev/null | wc -l)
		    if [ $npld -ne 1 ]; then
                    	echo "Found $npld plumed.dat files in dir `pwd`" >&2; exit 1
		    fi
		    plumed="-plumed $(ls *.dat)"   
		fi

                mdrunCall="$mpirun$mdrun"
                mdrunArgs="$cpiArg -stepout $stepout $verb -s $tpr $maxh $plumed $dd $npme $deffnm $opt $edi $pinArgs $gpuID_flag"
                echo -e "$mdrunCall $ntFlag $mdrunArgs >& md$key.lis &\n"
                let idir++
            done
            echo 'wait'
            cd $dir
        fi
        # echo "$scratchend"
    fi
} >> $dir/$jobfile

if [ $bGo = 1 ];then
    case $Qsystem in
        slurm)
            sbatch $jobfile >& sbatch${key}.out
            cat sbatch${key}.out >&2
            ;;
        *)
            echo "Unknown Qsystem = $Qsystem"; exit 1
            ;;
    esac
fi

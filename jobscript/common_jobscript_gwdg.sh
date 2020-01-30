#!/bin/bash
shopt -s -o nounset

=() { awk 'BEGIN {print '"$*"'}' < /dev/null; }
round_ndig () { awk '{printf "%.'$1'f\n", $1}'; }

version_check=""
version=''
d=0
m1=0
h=0
m2=0
jobname=cmbjob
tpr=""
npme=""
dd=""
deffnm=""
nnodes=1
ppn=not-set
md=1
exe=""
opt=""
mpi=0
gmxrc=gmxrc_not_set
mdrun=/usr/users/cmb/shared/bin/run_mdrun.sh
bGaussian=0
verb="-v"
bGo=0
key=""
email="${USER}@gwdg.de"
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
nt=""
nameext=""
ptile=""
np=""
ncores=""
bEmpty=0
bLatest=0
excludeNodes=''
bEnsureFullNode=0  # if you run 8 2-core jobs (with pinning, make sure the node is filled by your jobs): #BSUB -R np16                                                                            
batchInitLine=''

sbatch_tempfile=`mktemp sbatch.tempXXXXX`
#rm -f $sbatch_tempfile
trap "rm -rf $sbatch_tempfile" EXIT

loc=0
copyt=0

case `hostname` in
    gwdu102|gwdu103)
        ppn=16
        gmxrc=/home/uni05/cmb/shared/opt/gromacs/sandy-bridge/4.65/bin/GMXRC.bash
        queue=gpu-hub
        machine=sandybridge
        ;;
    smaug)
	ppn=12
	gmxrc=/data/shared/opt/gromacs/2018.6/bin/GMXRC.bash
	queue=gpu-hub
	machine=smaug
	;;
    *)
        echo -e "WARNING: Default ppn and gmxrc unknown (hostname = `hostname`).\n"
        exit 1
esac


#Set defaults
case $machine in
    sandybridge|broadwell)
        Qsystem=slurm
        walltime='2-00:00'
	intime=48
        maxh=48
	deplineTempl='#SBATCH --depend=afterok:JOBID'
	Qkey='#SBATCH'
        ;;
    *)
        echo "Unknown machine $machine"; exit 1
esac

echo "Found machine = $machine -- Qsystem $Qsystem, default walltime $walltime"


mdrun=run_mdrun.sh
mdrun_line=""
ntFlag=""
bNoNtFlag=0
bCpi=1
scratch=0
bReqScratch=0
gpuGeneration='pascal40'


while [ $# -gt 0 ]; do
    case "$1" in
	-empty) bEmpty=1;;
	-go) bGo=1;;
	-p) shift
            queue=$1
	    if [ $queue == 'medium' ]; then
		med='#SBATCH -A all'
	    else
		med=''
	    fi
	    ;;
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
	-loc)
	    loc=1
	    echo loc = $loc
	    copyt=1 ;;
	-copyt)
	    shift
	    copyt=$1 ;;
	--time|-t)
	    shift
	    intime=$1
	    d=$(= "int(($intime)/24)")
	    h=$(= "int(($intime)%24)")
	    m1=$(= "(($intime)%24)-$h")
	    m2=$(= "int($m1*60)")
	    if [ $Qsystem = slurm ]; then
		#slurm want D-HH:MM
		walltime="$d-$h:$m2"
	    else
		echo "Unknown Qsystem = $Qsystem"; exit 1
	    fi
	    ;;
	--jobname|-J|-j)
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
	    #echo $2
      	    #if [[ $2 == ^[0:9]* ]]; then
	    #  	shift
            #   maxnodes=$1
            #   echo maxnodes=$1
	    #else
            #   echo "No range of nodes given. Using $nnodes."
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
        -mdrun_line) shift
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
        -gpu-generation|-gpu-gener)
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
	-batch_line)
	    shift
	    batchInitLine=$1;;
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
        -scratch)
            scratch=1 ;;
        -req-scratch)
            bReqScratch=1 ;;
        -exclude-nodes)
            shift
            excludeNodes="$(echo "$1" | sed 's/@/ /g')" ;;
	-version)
            shift
	    version_check="1"
            version=$1
	    if [[ `hostname` == "gwdu103" ]]
	    then
		gmxrc="/usr/users/cmb/shared/opt/gromacs/broadwell/$version/bin/GMXRC"
	    elif [[ `hostname` == "gwdu102" ]]
	    then
		gmxrc="/usr/users/cmb/shared/opt/gromacs/sandy-bridge/$version/bin/GMXRC"
	    fi
            ;;
        *)
            echo -e "\n$0: Error, unknown argument: $1"
            exit 192
            ;;
    esac
    shift
done

echo verb="$verb"

[ $bCpi = 1 ] && cpiArg="-cpi" || cpiArg=""

if [ "$loc" -eq 1 ]; then
    [ $bCpi = 1 ] && cpiArg="-cpi ${dir}/state.cpt" || cpiArg=""
    maxh=$(echo "$intime-$copyt" | bc)
    if [ "$(echo "$maxh <= 0" | bc)" -eq 1 ]; then
        echo "you won't have enough time to copy files from the local disk, change the -copyt or increase -t option !!!"
        exit -1
    fi
fi
echo walltime=$walltime
echo maxh="$maxh"


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

if echo $mdrun | grep -q run_mdrun.sh; then
    gmxrcLine=""
else
    gmxrcLine="source gmxrc_and_modules.sh"
fi

if [[ $version_check = "1" ]]; then
    gmxrcLine="source $gmxrc"
fi


# edi present?
nedi=$(ls *edi 2>/dev/null | wc -l)
if [ $nedi -eq 1 ]; then
    edi="-ei $(ls *edi)"
elif [ $nedi -eq 0 ]; then
    edi=""
else
    echo "There are $nedi edi files in current dir. ??" >&2
fi

if [ $scratch = 1 ] ; then
    scratchbeg="                                                                                           
export WORK_FS=$dir                                                                                        
export JOB_FS=/scratch/${USER}/\$LSB_JOBID                                                                 
echo Job-Directory \$JOB_FS                                                                                
                                                                                                           
mkdir \$JOB_FS                                                                                             
cp ./* \$JOB_FS                                                                                            
cd \$JOB_FS                                                                                                
"
    scratchend="                                                                                           
cp ./* \$WORK_FS                                                                                           
cd \$WORK_FS                                                                                               
rm -rf \$JOB_FS                                                                                            
"
else
    scratchbeg=""
    scratchend=""
fi


ldPath=""
initMPI=""
spanline=""
#batchInitLine=""
if [ "$np" = "" ]; then
    np=$[nnodes*ppn]
fi


ldPath=''

case $machine in
    broadwell|sandybridge)                                                
        #[ $ptile = unset ] && ptile=$ppn
        #spanline="#SBATCH --tasks-per-node=[$ptile]"
        if [ $bMPI = 1 ]; then
            mpirun="mpirun  -np $np "
            ldPath="$ldPath:"
            mdrun="${mdrun}_mpi"
        else
            mpirun=""
        fi
        ;;
    interlagos)
        [ $ptile = unset ] && ptile=$ppn
        spanline="#SBATCH --tasks-per-node=[$ptile]"
        ldPath=""
        if [ $bMPI = 1 ]; then
            mpirun="mpirun.slurm -np $np "
	        # next line commented out, because apparently its not needed in slurm
            # echo "#SBATCH -a intelmpi" >>$sbatch_tempfile
            # ldPath="$ldPath:/cm/shared/apps/openmpi/gcc/64/1.6.4/lib64"                                                    
            ldPath=""
            mdrun="${mdrun}_mpi"
        else
            mpirun=""
        fi
        # make sure we pin the threads if running multiple mdrun on one node
        # However, if the node is not full, we should not pin. Otherwise, different                                                       
        # mdrun may use the same nodes (fixed below)                                                                                      
        ;;
    *)
        echo -e "\nError: Unknown machine $machine"; exit 1
        ;;
esac

initLD='export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:'"$ldPath"

nGPUsPerNode=0
logicalCoresPerPhysical=1
if [ $queue = gpu ] || [ $queue = gpu-hub ]; then
    if [ "$ncores" = "" ]; then
        ncores=$ppn
    fi

    case "$gpuGeneration" in
        kepler)
            # Use Kepler GTX 1070: nvgen=3D2    
            nGPUsPerNode=1
            logicalCoresPerPhysical=2
            {
                echo '#SBATCH --gres=gpu:gtx1070:1'                                                                                
		echo '#SBATCH --exclude=dge[015-045]'
                queue=gpu-hub                                                                                                                                                          
            } >> $sbatch_tempfile
            ;;
        maxwell)
            # Use GTX 980                                                                                                             
            nGPUsPerNode=2
            {
                echo '#SBATCH --gres=gpu:gtx980:1'
		echo '#SBATCH --exclude=gpu-hub' # Exclude the our own Pascal nodes
            } >> $sbatch_tempfile
            ;;
	maxwell4)
            # Use nodes with 4 GTX 980
            nGPUsPerNode=4
            {
                echo '#SBATCH --gres=gpu:gtx980:1'
            } >> $sbatch_tempfile
            ;;
        tesla)
            # Use the expensive Tesla (allows double prec.)
            nGPUsPerNode=2
            {
		#NOT IMPLEMENTED YET
                echo "tesla is not supported at the moment"
                
            } >> $sbatch_tempfile
            ;;
        gtx1080)
            # Use GTX 1080                                                                                                      
            nGPUsPerNode=2
            {
                echo '#SBATCH --gres=gpu:gtx1080:1'
		echo '#SBATCH --exclude=gpu-hub' # Exclude our own Pascal nodes
            } >> $sbatch_tempfile
            ;;
        pascal40)
            # Use GTX 1080                                                                                              
            nGPUsPerNode=4
            logicalCoresPerPhysical=2
            {
                echo '#SBATCH --gres=gpu:gtx1070:1'
                queue=gpu-hub
            } >> $sbatch_tempfile
            ;;
        *)
            echo "ERROR, unknown GPU generation: \"$gpuGeneration\". Uset option -gpu-gener or -gpu-generation" >&2; exit 1
    esac

    if [ "$ptile" = "" ]; then
        spanline="#SBATCH --ntasks-per-node=1"
    else
        spanline="#SBATCH --ntasks-per-node=$ptile"
    fi

    if [ $gpuGeneration = pascal40 ]; then
        # Jochen, Jun 16, 2017
        # Echo work-around while ptile option on Pascal40 nodes behaves strange
        spanline="#SBATCH -N 1"
    fi

    if [ "$gpu_shares" = "" ]; then
        # if GPU shares not given, use the same as number of physical CPU cores requested.
        gpu_shares=$[ncores]
    fi
    echo "GPU gener = $gpuGeneration, ncores = $ncores, gpu_shares = $gpu_shares"
    # echo "#SBATCH --licenses=[ngpus_shared=10]" >> $sbatch_tempfile
    # ppn=""
fi



if [ "$multidirs" = "" ] ;then
    jobfile=job${key}.sh
else
    jobfile=job.${jobname}${key}.sh
    # Check if the node is really full
    [ "$nnodes" -gt 1 ] && { echo -e "\nnnodes = $nnodes: Multidirs does not make sense." >&2; exit 1; }
    [ "$nt"      = "" ] && { echo -e "\nWith multidirs you must specify nt." >&2; exit 1; }
    ndir=$(echo $multidirs | wc -w | awk '{print $1}')
    npNeeded=$[nt*ndir/logicalCoresPerPhysical]
    echo "Multiple mdrun per node: ndir = $ndir, need $npNeeded phys. cores"

    # Note: np       = number of physical cores available/requested by command line
    #       npNeeded = number of physical cores needed

    if [ $npNeeded -gt $np ]; then
        echo -e "\nError, your job will need $npNeeded physical cores per node (nt=$nt, dir=$ndir), but you requested only $np physical cores."; exit 1
    fi

    if [ $npNeeded -lt $np ]; then
        echo -e "\nNOTE: Running $ndir jobs (nt=$nt) on $np cores would wastes recources. Reducing # of allocated cores and ptile to $npNeeded" >&2
        echo -e "NOTE: Will not pin mdrun threads to cores." >&2
        np=$npNeeded
        spanline="#BSUB -N $npNeeded"
        bPin=0
    fi
fi



if [ -e "$jobfile" ]; then
    echo Backing up $jobfile to $jobfile.bak
    mv $jobfile $jobfile.bak
fi



############################################################################################################################
                    
############################################################################################################################

[ "$bEmpty" = 1 ] && gmxrcLine=""
[ "$ncores" = "" ] && ncores=$np
[ "$bReqScratch" = 1 ] && echo "#SBATCH -C scratch" >> $sbatch_tempfile
[ "$spanline" != "" ]  && echo "$spanline"        >> $sbatch_tempfile
[ "$bLatest" = 1 ]     && echo "#SBATCH -C latest"  >> $sbatch_tempfile


if [[ "$Qsystem" = slurm ]]; then
    if [ $bEnsureFullNode = 1 ]; then
        if [[ ( ($ncores = 12) || ( $ncores = 16 ) || ( $ncores = 20 ) || ( $ncores = 24 ) ) && ( $queue = mpi ) ]]; then        
            echo "#SBATCH --batch=\"ncpus=${ncores}\"" >> $sbatch_tempfile
            echo "Will make sure that the node is full, using purely nodes with ${ncores} cores"
        fi
    fi
    {
	cat <<EOF
#!/bin/bash
#SBATCH -p $queue$qExtension
#SBATCH -o $dir/myjob${key}.out
#SBATCH -e $dir/myjob${key}.err
#SBATCH -c $(echo $logicalCoresPerPhysical*$ppn | bc)
#SBATCH -t $walltime
#SBATCH --job-name=$jobname$key
#SBATCH --mail-user=$email
$batchInitLine
$depline
$med

EOF

	cat $sbatch_tempfile

	cat <<EOF
echo Starting job at \$(date)
echo Jobid \$SLURM_JOBID
echo Host \$SLURM_JOB_NODELIST
echo Jobname \$SLURM_JOB_NAME
echo Subcwd \$SLURM_SUBMIT_DIR

cd $dir
EOF
    } > $jobfile

if [ $loc = 1 ]; then
    {
        echo -e "mkdir /local/${USER}_\$SLURM_JOB_ID\ncd /local/${USER}_\$SLURM_JOB_ID"
    } >> $jobfile
fi

    {
	cat <<EOF
$gmxrcLine
$line
$initLD
$scratchbeg

EOF
    } >> $jobfile

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


{

    if [ "$mdrun_line" != "" ]; then
        
        echo "Found complete mdrun command line (option -mdrun_line), do not write mdrun line automatically" >&2
        echo "$mdrun_line &> md${key}.lis"
        
    else
        
        pinArgs=""
        if [ "$multidirs" = "" ]; then
            [ "$bPin" = 1 ] && pinArgs="-pin on -pinoffset 0 -pinstride 1"
            if [[ "$md" = 1 ]]; then
                if [[ ( $queue = gpu ) && ( "$gpu_id" != unset ) ]]; then
                    [ "$gpu_id" = '' ] && gpuID_flag="-gpu_id 0" || gpuID_flag="-gpu_id $gpu_id"
                else
                    gpuID_flag=""
                fi
                
                mdrunCall="$mpirun$mdrun"
                mdrunArgs="$cpiArg -stepout $stepout $verb -s "${dir}/${tpr}" -maxh $maxh $dd $npme $deffnm $opt $edi $pinArgs"
		if [ "$loc" = 0 ]; then
                    echo "$mdrunCall $ntFlag $mdrunArgs $gpuID_flag >> md${key}.lis 2>&1"
		else
		    echo "$mdrunCall $ntFlag $mdrunArgs $gpuID_flag >> /local/${USER}_\$SLURM_JOB_ID/md${key}.lis 2>&1"
		fi
            else
                echo "$mpirun" "$exe"
            fi
        else
            idir=0
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
                    pinArgs="-pin on -pinoffset $[idir*nt] -pinstride 1"
                fi
                if [ $queue = gpu ]; then
                    # Get number of mdruns running per GPU. Round up, important if ndir is an odd number
                    nDirPerGPU=$(= "$ndir/$nGPUsPerNode+0.01" | round_ndig 0)
                    # nDirPerGPU=$[ndir/nGPUsPerNode]
                    gpuID_flag="-gpu_id $[idir/nDirPerGPU]"
                else
                    gpuID_flag=""
                fi
                
                mdrunCall="$mpirun$mdrun"
                mdrunArgs="$cpiArg -stepout $stepout $verb -s $tpr $maxh $dd $npme $deffnm $opt $edi $pinArgs $gpuID_flag"
                echo -e "$mdrunCall $ntFlag $mdrunArgs >& md${key}.lis &\n"
                let idir++
            done
            echo 'wait'
            cd $dir
        fi
        echo "$scratchend"
    fi
} >> $dir/$jobfile

if [ $loc = 1 ]; then
    {
        echo -e "echo \"copying files from /local/${USER}_\${SLURM_JOB_ID} to \${SLURM_SUBMIT_DIR} at \$(date)\""
	echo -e "cp --backup --suffix=.old_key${key} /local/\${USER}_\${SLURM_JOB_ID}/* \${SLURM_SUBMIT_DIR}"
	echo -e "echo \"Ending job at \$(date)\""
    } >> $jobfile
fi


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

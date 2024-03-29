#!/bin/bash
shopt -s -o nounset

sbatch2id(){
    tail -n1 | awk '{print $4}' | sed 's/[<,>]//g'
}

nChain=1


loc=1

if [ $# -gt 0 ]; then
    arg1=$1
else
    arg1=""
fi

if ! which jobscript_slurm.sh >& /dev/null; then
    echo jobscript.sh not in path.
    exit 1
else
    script=`which jobscript_slurm_gwdg.sh`
    echo -e "Using this script: $script\n"
fi

if [ "$arg1" != '-n' ]; then
    echo -e "This script calls `which jobscript.sh` multiple times to submit a job chain."
    echo -e 'Hint: First call jobscript.sh and check if the "job.sh" looks fine.\n'
    echo -e "Usage:\n$0 -n N [-go] [jobscript.sh-options]\n"
    exit 2
fi

log=chain.log
rm -f $log

# nr of jobs in chain
call="$0 $@"
nChain=$2
shift 2
echo -e "Generate chain of $nChain jobs.\n"
echo -e "Used this command:\n${call}\n" | tee -a $log

if [ "$1" = -ext ]; then
    log=chain${2}.log
    shift 2
else
    log=chain.log
fi

# arguments for jobscript.sh
opt="$@"

bGo=0
if [ $# -gt 0 ]; then
    echo $opt | grep -wq -- -go && bGo=1
    echo "Found bGo = $bGo"
fi

echo "Starting job chain of $nChain jobs at `date` by user $USER" >> $log
for ((i=0; i<nChain; i++)); do

    if [ $i == 0 ]; then
	echo "$script -key .ch$i $opt"
	$script -key .ch$i $opt || exit 1
    else
	echo "$script -key .ch$i -dep $idprev $opt"
	$script -key .ch$i -dep $idprev $opt || exit 1
    fi  
    
    {
    	cat <<EOF
# echo "some noncence every .sh"
EOF
        } >>*.ch$i.sh


    # pick jobID (use a name id_jobX for testing)
    if [ $bGo = 1 ]; then
	if [ -e sbatch.ch$i.out ]; then
	    idprev=$(sbatch2id < sbatch.ch$i.out)
	elif [ -e msub.ch$i.out ]; then
	    idprev=$(cat msub.ch$i.out | tail -n 1)
	else
	    echo -e "\nERROR: Did not find sbatch.ch$i.out or msub.ch$i.out in dir `pwd`\n" >&2; exit 1
	fi
	# Check if $idprev is really a number (and nothing more)
	if ! echo "$idprev" | egrep -q '^[0-9]+$|hannover.[0-9]+'; then
	    echo "Error, job ID ($idprev) seems incorrect."; exit 1
	fi
    else
	idprev="id_job$i"
    fi
    echo "Submitted job.ch$i.sh - jobID $idprev" | tee -a chain.log
done


# some things to put at the last member of the chain like concatenation of the files
#some if statement to check if local is active
{
	cat <<EOF
echo "last of the chain"
EOF
    } >>*.ch$(echo "$nChain-1"| bc ).sh


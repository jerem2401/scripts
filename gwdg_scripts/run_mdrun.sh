#!/bin/bash

S=$(basename $0)

# default version
gmxVersion=2016.3

# gromacs installation dirs
d0=/usr/users/cmb/shared/opt/gromacs

# mdrun suffix
suff=""

if [ "$1" = "-version" ]; then
    gmxVersion=$2
    shift 2
fi
if [ "$1" = "-mpi" ]; then
    suff=_mpi
    shift
fi

#arch=$(/home/uni05/cmb/shared/bin/gwdg-arch.sh) || exit 1
arch=$(/usr/users/jlapier/bin/gwdg-arch.sh) || exit 1

module load rev/20.12

if [ -e "$d0/$arch/$gmxVersion/gwdg_modules.list" ]; then
    cat $d0/$arch/$gmxVersion/gwdg_modules.list > req-modules.tmp
    echo "Modules: taken from $d0/$arch/$gmxVersion/gwdg_modules.list"
else
    # load required modules
    echo "Modules: defaults from gmx_required_modules.sh"
    #/home/uni05/cmb/shared/bin/gmx_required_modules.sh > req-modules.tmp
    /usr/users/jlapier/bin/gmx_required_modules.sh > req-modules.tmp
fi

for i in `cat req-modules.tmp`; do
    echo "$S: loading module $i"
    module load "$i"
done

# remove names behind the minus, such as 5.14-something to 5.14
gmxVersionNr=$(echo $gmxVersion | cut -d- -f1)

if [ -e $d0/$arch/$gmxVersion/bin/gmx${suff} ]; then
    executable=$d0/$arch/$gmxVersion/bin/gmx${suff}
    mdrun="$d0/$arch/$gmxVersion/bin/gmx${suff} mdrun"
elif [ -e $d0/$arch/$gmxVersion/bin/mdrun${suff} ]; then
    mdrun=$d0/$arch/$gmxVersion/bin/mdrun${suff}
    executable=$mdrun
else
    echo -e "run_mdrun.sh: ERROR, neither gmx${suff} not mdrun${suff} found in directory: $d0/$arch/$gmxVersion/bin/" >&2
    echo "Varialbes: gmxVersion = $gmxVersion -- arch= $arch -- suff = $suff" >&2
    exit 1
fi

ncores=$(grep -c "^processor" /proc/cpuinfo)
echo "$S: Node = "$(hostname)
echo "$S: Architecture = $arch"
echo "$S: $ncores CPU-cores"
echo "$S: mdrun = $mdrun"
echo

$mdrun $@

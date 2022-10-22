
version=4.68

if [ $# -gt 0 ]; then
    if [ "$1" = -version ]; then
	version=$2
    else
	version=4.68
    fi
fi

#gwdg_arch=$(/home/uni05/cmb/shared/bin/gwdg-arch.sh)
gwdg_arch=$(/usr/users/jlapier/bin/gwdg-arch.sh)
echo "gmxrc_and_modules.sh: Sourcing Gromacs, architecture $gwdg_arch, version $version"

#source modules
if [ -e "/usr/users/cmb/shared/opt/gromacs/$gwdg_arch/$version/gwdg_modules.list" ]; then
    modules=$(cat /usr/users/cmb/shared/opt/gromacs/$gwdg_arch/$version/gwdg_modules.list)
    echo "Modules taken from /usr/users/cmb/shared/opt/gromacs/$gwdg_arch/$version/gwdg_modules.list"
else
    # load default modules
    echo "Modules: defaults from gmx_required_modules.sh"
    #modules=$(/home/uni05/cmb/shared/bin/gmx_required_modules.sh)
    modules=$(/usr/users/jlapier/bin/gmx_required_modules.sh)
fi

echo "gmxrc_and_modules.sh: loading modules: $modules"
module load $modules


# source gromacs

load_gmxrc=/usr/users/cmb/shared/opt/gromacs/$gwdg_arch/$version/bin/GMXRC.bash
if [ ! -e "$load_gmxrc" ]; then
    echo -e "\ngmxrc_and_modules.sh: ERROR, not found: $load_gmxrc\n" >&2
fi
source $load_gmxrc

which_mdrun=$(which mdrun 2>/dev/null)
if [ "$which_mdrun" = "" ]; then
    which_mdrun=$(which gmx 2>/dev/null)
fi
echo "gmxrc_and_modules.sh: now mdrun in path: $which_mdrun"

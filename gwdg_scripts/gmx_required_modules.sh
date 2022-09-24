#!/bin/bash

#arch=$(/home/uni05/cmb/shared/bin/gwdg-arch.sh) || exit 1
arch=$( /usr/users/jlapier/bin/gwdg-arch.sh) || exit 1
case $arch in
    nehalem|gtx480|gpu|gtx770)
	# echo "intel/mkl/64/11.2/2015.3.187 intel-mpi/64/5.0.3/048 intel/compiler/64/15.0/2015.3.187 cuda75/toolkit/7.5.18"
	# echo "intel/mkl/64/11.2/2015.5.223 intel/compiler/64/15.0/2015.5.223 intel-mpi/64/5.1.2/150 cuda80/toolkit/8.0.44"

	# Since June 2017:
	# echo "intel/mkl/64/11.3.3/2016.3.210 intel/compiler/64/16.0.3/2016.3.210 intel-mpi/64/5.1.2/150 cuda80/toolkit/8.0.44"
	# Since July 2017:
	echo "intel/mkl/64/2017/2.174 intel/compiler/64/16.0.4/2016.4.258 intel/mpi/64/2017/2.174 cuda10.2/toolkit/10.2.89"
	;;
    broadwell|sandy-bridge)
	#No open mpi because openmpi/3.1.4 not compatible with gcc8.4.0 (to double check) and gcc/9.3.0 not compatible with
	#cuda10.2.89 but with cuda 11.2.0 and not sure if we want cuda 11 here for all versions.
	echo "cuda/10.2.89 gcc/8.4.0"
	;;
    *)
	echo "$(basename $0): invalid architecture: $arch" >&2
	;;
esac

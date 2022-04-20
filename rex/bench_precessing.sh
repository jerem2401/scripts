#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes


module load rev/20.12
module load anaconda3/2020.07
set +eu
source activate env1

sed -n -e '/avg:/,/std:/ p' rex_opti.out  | egrep -o '([0-9]*\.[0-9]*)' > avg.txt
sed -n -e '/std:/,/]/ p' rex_opti.out | egrep -o '([0-9]*\.[0-9]*)' > std.txt

bench_precessing.py > process_out.txt

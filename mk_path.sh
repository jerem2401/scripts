#!/bin/bash
###############################################################################
#Author: Jeremy Lapierre <jeremy.lapierre@uni-saarland.de>
#Option description:
#
###############################################################################
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # dont hide errors within pipes

pym=0
itraj="../nopbc2.xtc"
if [[ $HOSTNAME == gwdu* ]]; then
    ndx="/usr/users/jlapier/simulation/syncsim/pol/heavy_h/transfer/index_notwholedna.ndx"
    rndx="/usr/users/jlapier/simulation/syncsim/pol/heavy_h/ref/idx_notwhole.txt"
else
    ndx="/data/users/jeremy/simulation/syncsim/pol/heavy_h/ref/index_notwholedna.ndx"
    rndx="/data/users/jeremy/simulation/syncsim/pol/heavy_h/ref/idx_notwhole.txt"
fi
check_path=''

while [ $# -gt 0 ]; do
    case "$1" in
        -pym) pym=1;;
	-check) shift
	    check_path=$1;;
	-itraj) shift
	    itraj=$1;;
	-ndx) shift
	    ndx=$1;;
	-rndx) shift
	    rndx=$1;;
	-lmd) shift
	    lmd=$1;;
	*) echo -e "\n$0: Error, unknown argument: $1"
	   exit 192;;
    esac
    shift
done

if [ -z "$check_path" ]; then
    if [ $pym = 0 ]; then
        gmx trjconv -nice 0 -f $itraj -o testnotrn.pdb -n $ndx -s ../md.tpr
        plumed pdbrenumber --ipdb testnotrn.pdb --opdb test.pdb --atomnumbers $rndx
        sed '/TER/d;/MODEL/d;/^END$/d;/REMARK/d;/TITLE/d;/CRYST/d' test.pdb > ref_clean.pdb

        gawk '{
    pdb[NR]=$0
    split(pdb[FNR],flds,FS,seps)
        for (i=1;i in flds;i++){
            if (flds[3]!="CA"){
                flds[9]="0.00"
                flds[10]="1.00"
            }
        printf "%s%s", flds[i], seps[i]
        }
    print ""
    }' ref_clean.pdb > ref_clean_2.pdb
    else
        sed '/TER/d;/MODEL/d;/^END$/d' test.pdb > ref_clean.pdb 
        gawk '{
    pdb[NR]=$0
    split(pdb[FNR],flds,FS,seps)
        for (i=1;i in flds;i++){
            if (flds[3]!="CA"){
                flds[10]="0.00"
                flds[11]="1.00"
            }
        printf "%s%s", flds[i], seps[i]
        }
    print ""
    }' ref_clean.pdb > ref_clean_2.pdb
    fi

    csplit --suppress-matched ref_clean_2.pdb /ENDMDL/ '{*}' -b %02d.pdb
    wait
    last=$(command ls xx* -v | tail -1)
    rm ./$last
    wait

    echo "removing previous plumed.dat"
    rm -f ./plumed.dat

    for i in $(command ls -v xx*); do
        echo "RMSD REFERENCE=$i TYPE=OPTIMAL LABEL=RMSD${i:2:-4} NOPBC" >> plumed.dat
    done
    echo "PRINT ARG=* STRIDE=1 FILE=colvar_RMSDs.txt" >> plumed.dat
    plumed driver --plumed plumed.dat --mf_xtc ../nopbc2.xtc

    awk '(NR < 7) {print $1" "$2" "$3" "$4" "$5" "$6}' colvar_RMSDs.txt
else
    bck=$(command ls path_${check_path%.txt}.pdb &>/dev/null || unset)
    nbck=$(echo $bck | wc -w)

    if [ ! -z "$bck" ]; then
	path="path_${check_path%.txt}_${nbck}.pdb"
    else
	path="path_${check_path%.txt}.pdb"
    fi

    while read line; do
        for i in $line; do
	    cat "xx${i#RMSD}.pdb" >> $path
	    echo "END" >> $path
	done
    done <"$check_path"

    echo -e "p: PATHMSD REFERENCE=$path LAMBDA=$lmd  NOPBC\nPRINT ARG=* STRIDE=1 FILE=colv_${path%.pdb}${lmd}.txt" >> plumed_${path%.pdb}${lmd}.dat
    plumed driver --plumed "plumed_${path%.pdb}${lmd}.dat" --mf_xtc ../nopbc2.xtc
fi

#sed '/TER/d;/MODEL/d;/^END$/d' test.pdb > ref_clean.pdb
##############!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#################
#gawk -f ~/simulation/syncsim/pol/heavy_h/pull_phase3_ZN/my.awk ref_clean.pdb > ref_clean_2.pdb
#csplit --suppress-matched ref_clean_2.pdb /ENDMDL/ '{*}' -b %02d.pdb
#wait
#last=$(command ls xx* -v | tail -1)
#rm ./$last
#wait
###echo "WHOLEMOLECULES ENTITY0=86293-87474 ENTITY1=12747-13311 ENTITY2=38673-38896 ENTITY3=8127-8278 ENTITY4=10583-11589 ENTITY5=21788-21985 ENTITY6=87475-90408" > plumed.dat
#for i in $(command ls -v xx*); do echo "RMSD REFERENCE=$i TYPE=OPTIMAL LABEL=RMSD${i:2:-4} NOPBC" >> plumed.dat; done
#echo "PRINT ARG=* STRIDE=1 FILE=colvar_RMSDs.txt" >> plumed.dat
#plumed driver --plumed plumed.dat --mf_xtc ../nopbc2.xtc



#plumed driver --plumed plumed.dat --mf_pdb test.pdb


#for i in RMSD01 RMSD03 RMSD43 RMSD68 RMSD74 RMSD80 RMSD107 RMSD123 RMSD127 RMSD134 RMSD136 RMSD138 RMSD144 RMSD149 RMSD152 RMSD157 RMSD188 RMSD213 RMSD217 RMSD233 RMSD252 RMSD335 RMSD346 RMSD384 RMSD405 RMSD456 RMSD483 RMSD545 RMSD574 RMSD633 RMSD663 RMSD692 RMSD708 RMSD718 RMSD795 RMSD801 RMSD818 RMSD823 RMSD843 RMSD854 RMSD871 RMSD876 RMSD947 RMSD987 RMSD1018 RMSD1067 RMSD1081 RMSD1133 RMSD1162 RMSD1211 RMSD1256 RMSD1262 RMSD1281 RMSD1304 RMSD1340 RMSD1374 RMSD1393 RMSD1463 RMSD1488 RMSD1511 RMSD1601 RMSD1628 RMSD1750 RMSD1765 RMSD1769 RMSD1798 RMSD1876 RMSD2000; do cat "xx${i#RMSD}.pdb"; echo "END"; done >> path.pdb

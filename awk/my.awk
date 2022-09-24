#replace a field (here field 10) from one pdb to entry in another one while keeping formatting.
#usage: gawk -f my.awk t1.txt t2.txt > out
#NR is total record number (total lines read), FNR is the record number (line number)
#If actual line number = total number of line then execute the following block
NR==FNR { pdb[NR]=$0; next }
{
split(pdb[FNR],flds,FS,seps)
flds[10]=$10
for (i=1;i in flds;i++)
printf "%s%s", flds[i], seps[i]
print ""
}

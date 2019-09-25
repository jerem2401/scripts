for i in rep*
do
	cd $i
	cpumb_plot.py -f eq2ax/colvar_pull.txt ax2eq/colvar_pull.txt --a 0.5
	echo "${i}-done"
	cd ..
done

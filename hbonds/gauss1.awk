# Script from Jochen S. Hub#
############################
#usage:
#head shr_wat-dna.csv -n1
#>1.9376  -46820.0        150.0
#awk -v sig=3 -f gauss1.awk < shr_wat-dna.csv > shr_wat-dna_smooth.dat
#head shr_wat-dna_smooth.dat
#>1.9376 -46506.3
############################
# simple Gauss filter for 1D data
# 
# expects xy input with SORTED x
#
# sig : width (sigma) of the filter (in units of x)
# win : maximal # of datapoints  to take into account on both sides
#       (everything within +/- ~3*sigma sould do)
#
BEGIN {win=60; two_sig2=2*sig*sig;};
{
    # Read in data
    x[NR-1]=$1
    y[NR-1]=$2
}
END{
    # maybe certain x-poisitions appear multiple times, so find first x[i] that differs from x[0]
    # normally, dx is simply x[1]-x[0], of course
    j = 1
    do{
	dx           = (x[j]-x[0])/j
	j++
    } while (dx == 0)
    threeSigBins = int(3*sig/dx)
    # print dx, threeSigBins
    
    for (i=0; i<NR; i++)
    {
        i0 = i-threeSigBins
        i1 = i+threeSigBins

        if (i0< 0)     i0 = 0
        if (i1>(NR-1)) i1 = NR-1

        W   = 0
        sum = 0
	cen = x[i]
        for (j = i0; j<=i1; j++)
        {
	    exp_arg = -(x[j]-cen)*(x[j]-cen)/(two_sig2)
	    if (exp_arg < -100)
		w = 0
	    else
		w = exp(exp_arg)
            sum += w*y[j]
            W   += w
        }

	# If the same x appears multiple times, only write the first occurance
	# For normal data sets, this has no effect
	if (i == 0 || x[i] != lastx)
	{
	    print x[i], sum/W
	}
	lastx = x[i]
    }
}

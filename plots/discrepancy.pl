# To see this plot, start gnuplot and type:
#
#   load "discrepancy.pl"

set autoscale
set logscale y

set title "L2 star discrepancy"
set xlabel "Number of points"
set ylabel "L2 star discrepancy"

#plot "../data/discrepancies/output.dat" every 10 with lines
plot "$0" with lines, "$1" with lines

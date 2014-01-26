# To see this plot, start gnuplot and type:
#   call "rel_err.pl" "../data/simulations/mc_result.dat"

set logscale y
set logscale x
set title "Relative error"
set xlabel "Number of points"
set ylabel "Abs(100*(Res-Exact)/Exact"
set label "dim=100" at graph 0.05,0.25
set label "Slope=0.21" at graph 0.05,0.2
set label "u(i) all different pseudo-random numbers between 0 and 1" at graph 0.05,0.15
set label "square root started from n=50002, and skipped the first 500 primes" at graph 0.05,0.1
set label "NbRuns=10" at graph 0.05,0.05
plot "../data/simulations/pseudorandom.dat" using 1:4 title "Pseudorandom" with lines , \
     "../data/simulations/square_root.dat" using 1:4 title "Square Root" with lines, \
     "../data/simulations/sobol.dat" using 1:4 title "Sobol" with lines, \
     1/sqrt(x), \
     1/x
set terminal postscript eps color
set out 'relative_error.eps'
replot
set term x11
#set terminal postscript eps color lw 15 "Helvetica" 20

# TODO: also do fitting and plot -alpha parameter!

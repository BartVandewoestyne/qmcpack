# To see this plot, start gnuplot and type:
#   call "abs_err.pl" "../data/simulations/mc_result.dat"

set logscale y
set logscale x
set title "Absolute error"
set xlabel "Number of points"
set ylabel "Absolute error"
plot "../data/simulations/pseudorandom.dat" using 1:3 title "Pseudorandom" with lines , \
     "../data/simulations/square_root.dat" using 1:3 title "Square Root" with lines, \
     "../data/simulations/sobol.dat" using 1:3 title "Sobol" with lines, \
     1/x
set key left bottom

# TODO: also do fitting and plot -alpha parameter!

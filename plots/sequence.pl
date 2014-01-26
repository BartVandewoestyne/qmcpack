# Plot a two-dimensional projection of the specified point set.

set terminal postscript eps enhanced
set out "myfile.eps"

set xlabel "Dimension X"
set ylabel "Dimension Y"
set title "MyTitle"
set xrange [0:1]
set yrange [0:1]
set size square
plot "myfile.dat" using 1:2 title "" with points pointtype 7

# To plot only lines 514 to 1025 of your data file:
#plot "$0" every ::513::1024 using $1:$2

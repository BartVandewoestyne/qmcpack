# Plot a two-dimensional projection of the specified point set.
# Plot point 1 up to point ...

set terminal postscript eps enhanced
set out "myfile.eps"

set xlabel "Dimension 1"
set ylabel "Dimension 2"
set title "MyTitle"
set xrange [0:1]
set yrange [0:1]
set size square
plot "myfile.dat" every ::0::endpoint-1 using 1:2 title "" with points pointtype 7
#     "myfile.dat" every ::endpoint-1::endpoint-1 using 1:2 title "" with points pointtype 7 ps 2, \

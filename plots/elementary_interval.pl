set terminal postscript eps enhanced
set out "myfile.eps"

b = 2
kx = 2
ky = 3

set title "MyTitle"
set xlabel "Dimension X"
set ylabel "Dimension Y"
set title "Elementary interval"
set xrange [0:1]
set yrange [0:1]
set size square
set xtics 1.0/b**kx
set ytics 1.0/b**ky
set grid front lt -1 lw 3
plot "myfile.dat" using 1:2 title "" with points pointtype 7

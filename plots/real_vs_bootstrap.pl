# Plot real data versus bootstrap estimated data for some process.


#set term postscript enhanced
#set out "myoutput.eps"

set xlabel "Time"
set ylabel "Location"

plot "mydata.dat" u 1:2 title "x_{real}" pt 6 w l, \
     "mydata.dat" u 1:4 title "x_{bootstrap}" pt 2 w l

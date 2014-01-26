#set term postscript enhanced
#set out "myoutput.eps"

set title "Discrete PDF at step mytimestep"
set xlabel "Samples"
set ylabel "Weights"

plot "mydata.dat" u samplecolumn:weightcolumn smooth frequency w p title ""
#plot "mydata.dat" u samplecolumn:weightcolumn title ""

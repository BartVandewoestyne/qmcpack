# Plot the value of the integral as the number of samples increases.


set autoscale

set xlabel "Number of points"
set ylabel "Integral value"

plot "myfile.dat" using mycolumn_n:mycolumn_res with lines

# This script plots a histogram of some given datafile, using a certain
# width for the bins.


# The width of the bins.
width = my_width

# A function that returns the nearest integer to a given number.
nint(x) = floor(x+0.5)

# The function transforming a point x to it's bin.
bin(x) = width*nint(x/width)

# Always start plotting y-values from zero.
set yrange [0:]

# Plot it!
plot "myfile.dat" using (bin($1)):(1.0) smooth frequency with histeps title ""

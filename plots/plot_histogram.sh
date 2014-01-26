#!/bin/bash

filename=$1     # e.g. file.dat
width=$2        # e.g. 0.2

# Note that we use another separator for the sed command, to eliminate the
# problem of having to escape the backslashes in the filename.
sed "s#my_width#$width#g; \
     s#myfile.dat#$filename#g" histogram.pl | gnuplot -persist

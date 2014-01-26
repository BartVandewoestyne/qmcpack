#!/bin/bash

filename=$1     # e.g. file.dat
column_n=$2     # e.g. 2 (the N-values are in the 2nd column)
column_res=$3   # e.g. 3 (the integral values for each N are in the 3rd column)

# Note that we use another separator for the sed command, to eliminate the
# problem of having to escape the backslashes in the filename.
sed "s#mycolumn_n#$column_n#g; \
     s#mycolumn_res#$column_res#g; \
     s#myfile.dat#$filename#g" integration.pl | gnuplot -persist

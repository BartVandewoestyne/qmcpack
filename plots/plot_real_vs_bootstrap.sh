#!/bin/bash

filename=$1             # e.g. file.dat
column_real=$2          # e.g. 2 (the real values are in the 2nd column)
column_bootstrap=$3     # e.g. 4 (the bootstrap estimated values for
                        #         each N are in the 4th column)

# Note that we use another separator for the sed command, to eliminate the
# problem of having to escape the backslashes in the filename.
sed "s#myoutput.eps#${filename//.dat/.eps}#g; \
     s#1:2#1:$column_real#g; \
     s#1:4#1:$column_bootstrap#g; \
     s#mydata.dat#$filename#g" real_vs_bootstrap.pl | gnuplot -persist

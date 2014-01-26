#!/bin/bash

filename=$1     # e.g. file.dat
timestep=$2     # e.g. 2 (we want the PDF at step 2)
nbsteps=$3      # e.g. 50 (there are 50 steps written to the file)

# Note that we use another separator for the sed command, to eliminate the
# problem of having to escape the backslashes in the filename.
sed "s#myoutput.eps#${filename//.dat/.eps}#g; \
     s#mytimestep#${timestep}#g; \
     s#samplecolumn#${timestep}#g; \
     s#weightcolumn#$((${timestep}+${nbsteps}))#g; \
     s#mydata.dat#$filename#g" discrete_pdf.pl | gnuplot -persist

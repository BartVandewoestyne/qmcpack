#!/bin/bash
#
# Generate incremental plots of the sequence.  That is: start plotting point 1,
# then add point 2, etc...

filename=$1     # e.g. ../tests/faure.dat
dim_x=$2        # e.g. 1
dim_y=$3        # e.g. 2
endpoint=$4     # e.g. 128
title=$5        # e.g. Halton

if [ "$#" != 5 ]; then

  echo -e "Usage: generate_incremental_sequence_plot.sh <file> <x_dimension> <y_dimension> <endpoint> <title>"
  echo -e "Example: generate_incremental_sequence_plot.sh test.dat 1 2 128 MyTitle"
  exit 1

else

for i in `seq 1 ${endpoint}`;
  do
    sed -e "s#MyTitle#${title}#g" \
        -e "s#myfile.dat#${filename}#g" \
        -e "s#myfile.eps#${filename%.dat}.eps#g" \
        -e "s#Dimension 1#Dimension ${dim_x}#g" \
        -e "s#Dimension 2#Dimension ${dim_y}#g" \
        -e "s#endpoint#${i}#g" \
        -e "s#1:2#${dim_x}:${dim_y}#g" incremental_sequence.pl | gnuplot
    gv ${filename%.dat}.eps
  done

fi

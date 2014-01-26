#!/bin/bash

filename=$1     # e.g. ../tests/faure.dat
dim_x=$2        # e.g. 1
dim_y=$3        # e.g. 2
title=$4        # e.g. Halton

if [ "$#" != 4 ]; then

  echo -e "Usage: generate_sequence_plot.sh <file> <x_dimension> <y_dimension> <title>"
  echo -e "Example: generate_sequence_plot.sh test.dat 1 2 MyTitle"
  exit 1

else

  sed -e "s#MyTitle#${title}#g" \
      -e "s#myfile.dat#${filename}#g" \
      -e "s#myfile.eps#${filename%.dat}.eps#g" \
      -e "s#Dimension X#Dimension ${dim_x}#g" \
      -e "s#Dimension Y#Dimension ${dim_y}#g" \
      -e "s#1:2#${dim_x}:${dim_y}#g" sequence.pl | gnuplot

fi

# TODO
#!/bin/bash

filename=$1     # e.g. ../tests/faure.dat
dim_x=$2        # e.g. 1
dim_y=$3        # e.g. 2
base=$4         # The base of the elementary interval
d_x=$5          # The x-power (fineness) of the elementary interval (>=0)
d_y=$6          # The y-power (fineness) of the elementary interval (>=0)
title=$7        # e.g. Halton

if [ "$#" != 7 ]; then

  echo -e "Usage: generate_elementary_interval_plot.sh <file> <x_dimension> <y_dimension> <base> <d_x> <d_y> <title>"
  echo -e "Example: generate_elementary_interval_plot.sh test.dat 1 2 2 1 2 MyTitle"
  exit 1

else

  sed -e "s#MyTitle#${title}#g" \
      -e "s#myfile.dat#${filename}#g" \
      -e "s#myfile.eps#${filename%.dat}.eps#g" \
      -e "s#Dimension X#Dimension ${dim_x}#g" \
      -e "s#Dimension Y#Dimension ${dim_y}#g" \
      -e "s#b = 2#b = ${base}#g" \
      -e "s#kx = 2#kx = ${d_x}#g" \
      -e "s#ky = 3#ky = ${d_y}#g" \
      -e "s#1:2#${dim_x}:${dim_y}#g" elementary_interval.pl | gnuplot

fi

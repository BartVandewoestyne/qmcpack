#!/bin/bash

filename=$1     # e.g. ../tests/faure.dat
max_dim=$2

if [ "$#" != 2 ]; then

  echo -e "Usage: show_all_projections.sh <file> <max_dim>"
  echo -e "Example: show_all_projections.sh test.dat 10"
  exit 1

else

dim_x=1

while [[ $dim_x -le $max_dim && "$answer" != "n" ]]; do

  dim_y=$(($dim_x+1))
  while [[ $dim_y -le $max_dim && "$answer" != "n" ]]; do

    sed -e "s#MyTitle##g" \
        -e "s#myfile.dat#${filename}#g" \
        -e "s#myfile.eps#${filename%.dat}.eps#g" \
        -e "s#Dimension X#Dimension ${dim_x}#g" \
        -e "s#Dimension Y#Dimension ${dim_y}#g" \
        -e "s#1:2#${dim_x}:${dim_y}#g" sequence.pl | gnuplot

    echo "X = $dim_x  Y = $dim_y"

    if [[ $dim_x == 1 && $dim_y == 2 ]]; then
      gv --watch ${filename%.dat}.eps &
    fi

    echo -n "Continue (y/n)? "
    read answer

    dim_y=$(($dim_y+1))

  done

  dim_x=$(($dim_x+1))

done

rm ${filename%.dat}.eps

fi

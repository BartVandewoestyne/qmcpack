#!/bin/bash
# 2008-11-06 [koenp] script added to cvs-tree

SRC_DIR=../../src

## update dependencies
pushd .
cd $SRC_DIR
make depend.inc
popd

## parse them to dot
./dependtodot.sh "mod_cli,mod_bootstrap_filter,mod_system_model,mod_measurement_model" < "$SRC_DIR/depend.inc" > depend_full.dot

## remove transitive dependencies
cp depend_full.dot depend.dot
cp depend.dot depend.tmp ; tred depend.tmp > depend.dot
cp depend.dot depend.tmp ; unflatten -l 2 -f -c 2  depend.tmp > depend.dot

## Create figure
dot -Tpdf -odepend.pdf depend.dot 

## clean up

rm -f depend_full.dot depend.tmp depend.dot

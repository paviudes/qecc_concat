#!/bin/bash
make clean
make
echo "Moving bmark.so to chflow/src/simulate"
mv bmark.so ./../chflow/src/simulate/
ls -lt ./../chflow/src/simulate/bmark.so

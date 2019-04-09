#!/bin/bash
#PBS -l walltime=1000000000

cd /exports/pierre/cbird

python2 mcbird.py $aa $b $c $d $e $f

wait

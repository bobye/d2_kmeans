#!/bin/sh

#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -l pmem=1gb
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8


## run the 0-th phase on 2000 image data with 10 clusters
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 8,8 --phase_only 0 --clusters 10

## run the 1-th phase on 2000 image data with 10 clusters
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 8,8 --phase_only 0 --clusters 10

## run all phases on 5000 image with 10 clusters
time ./d2 -i data/icip14_data/total.d2 -p 2 -n 5000 -d 3,3 -s 8,8  --clusters 10



cd data/protein_seq/

## run 1-gram on protein data with 3 clusters
time ./protein 0 3

## run 2-gram on protein data with 3 clusters
time ./protein 1 3

## run 3-gram on protein data with 3 clusters
time ./protein 2 3


## run 1,2,3-gram on protein data with 2 clusters
time ./protein -1 2





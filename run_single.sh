#!/bin/sh

#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -l pmem=1gb
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8
time ./d2 -i mountaindat.txt -p 2 -n 2000 -d 3,3 -s 6,11 --phase_only 0 --clusters 10 --centroid_method 0

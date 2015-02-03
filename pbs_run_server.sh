#!/bin/sh

#PBS -l nodes=1:ppn=8
#PBS -l walltime=1:00:00
#PBS -l pmem=1gb
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=8

module load openmpi/gnu

time mpirun -n 5 ./d2 -i data/icip14_data/image.d2 -p 2 -n 1000 -d 3,3 -s 8,8 --phase_only 0 --clusters 10 -o centroids.d2 --centroid_method 0  --max_iter 20


num_of_nodes=4
batch_size=`mpirun -n 1 ./d2 -i data/icip14_data/total.d2 -p 2 -n 5000 -d 3,3 -s 8,8 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 10 -o centroids.d2 --centroid_method 0  --max_iter 20


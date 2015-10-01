#!/bin/bash

#PBS -l nodes=1:nogpu:ppn=32
#PBS -l walltime=1000:00:00
#PBS -l pmem=32gb
cd $PBS_O_WORKDIR

module load gcc-4.9.2
module load openmpi-1.6/gnu
main_dir=${PWD}

cd $main_dir
dataname=20news
num_of_nodes=32
s=${ss} # 32 64 100
dim=300
data_files=../docs_clustering_benchmark/$dataname\_cluster.d2s

# make sure split data only once
#batch_size=`./d2 -i $data_files -n 20000 -d $dim -s $s --types 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
batch_size=600

# run with different cluster numbers, use nohup to forward stdout
# results will save to different files with multiple runs
k=${kk} # 20 30 40 60 80
nohup mpirun -n $num_of_nodes ./d2 -i $data_files -n $batch_size -d $dim -s $s --clusters $k --types 7 > nohup.$dataname.$s.$k.out 

# save nohup.*.out for convergence analysis

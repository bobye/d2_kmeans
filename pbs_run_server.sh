#!/bin/sh

#PBS -l nodes=2:ppn=8
#PBS -l walltime=96:00:00
#PBS -l pmem=16gb
cd $PBS_O_WORKDIR

module load openmpi/gnu

#num_of_nodes=4
#batch_size=`./d2 -i data/icip14_data/total.d2 -p 2 -n 5000 -d 3,3 -s 8,8 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
#time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 10 -o centroids.d2 --centroid_method 0  --max_iter 20

num_of_nodes=16
data_files=data/20news_bydate/20newsgroups_clean/20newsgroups.d2s
k=10
batch_size=`./d2 -i $data_files -n 20000 -d 300 -s 128 --type 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
nohup mpirun -n $num_of_nodes ./d2 -i $data_files -n $batch_size -d 300 -s 64 --type 7 --clusters $k -o centroids.d2 > nohup.$num_of_nodes-$k.out


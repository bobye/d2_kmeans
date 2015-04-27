#!/bin/bash

#PBS -l nodes=1:ppn=8
#PBS -l walltime=96:00:00
#PBS -l pmem=8gb

cd $PBS_O_WORKDIR
module load openmpi/gnu
main_dir=${PWD}

ratios=(100 80 60 40)
ss=(80 64 48 32)

for i in 0 1 2 3; do
for k in 120 240 360; do
ratio=${ratios[$i]}
# extract data
cd data/usps
bunzip2 -k usps_blankout"$ratio"_train.d2s.bz2
bunzip2 -k usps_blankout"$ratio"_test.d2s.bz2
ln -s usps.d2s.vocab0 usps_blankout"$ratio"_train.d2s.vocab0
ln -s usps.d2s.vocab0 usps_blankout"$ratio"_test.d2s.vocab0

cd $main_dir
num_of_nodes=8
s=${ss[$i]}
train_files=data/usps/usps_blankout"$ratio"_train.d2s
test_files=data/usps/usps_blankout"$ratio"_test.d2s

# make sure split data only once
train_batch=`./d2 -i $train_files -n 8000 -d 2 -s $s --types 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
test_batch=`./d2 -i $test_files -n 3000 -d 2 -s $s --types 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`

# run with different cluster numbers, use nohup to forward stdout
# results will save to different files with multiple runs
nohup mpirun -n $num_of_nodes ./d2 -i $train_files -n $train_batch -d 2 -s $s --types 7 --clusters $k > nohup.$k.$ratio.out 

centroids_file=`grep "Write centroids to" nohup.$k.$ratio.out | head -1 | awk '{ print $NF }'`
	
# assign samples from test set to clusters
nohup mpirun -n $num_of_nodes ./d2 -i $test_files -n $test_batch -d 2 -s $s --types 7 --clusters $k --eval $centroids_file >> nohup.$k.$ratio.out

# save nohup.*.out for convergence analysis

done
done

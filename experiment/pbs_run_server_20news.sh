#!/bin/sh

#PBS -l nodes=2:ppn=8
#PBS -l walltime=96:00:00
#PBS -l pmem=16gb
cd $PBS_O_WORKDIR

module load openmpi/gnu
main_dir=${PWD}

# download data
cd data && git clone -o 20news-bydate https://github.com/bobye/20newsgroups.git 
cd 20news-bydate/20newsgroups_clean
bunzip2 20newsgroups.d2s.bz2
bunzip2 20newsgroups.d2s.vocab0.bz2


# run d2 clustering of documents
# make sure you change the line in src/d2/clusetering_io.c from 
#>   if (1) { // non pre-processed
# to 
#>   if (0) { // non pre-processed
# this will pre-process each BoW into a representation which has less than s supports.
# Of course, one can opt for pre-processing them on their own, but our implementation
# provides such function. 

cd $main_dir
num_of_nodes=16
s=16 # 2 4 8 16 32 64
data_files=data/20news-bydate/20newsgroups_clean/20newsgroups.d2s

# make sure split data only once
batch_size=`./d2 -i $data_files -n 20000 -d 300 -s $s --types 7 --prepare_batches $num_of_nodes --pre_process | grep batch_size | sed 's/^.*batch_size://g'`


# run with different cluster numbers, use nohup to forward stdout
# results will save to different files with multiple runs
k=10 # 20 30 40 60 80
nohup mpirun -n $num_of_nodes ./d2 -i $data_files -n $batch_size -d 300 -s $s --clusters $k > nohup.$k.out 

# save nohup.*.out for convergence analysis

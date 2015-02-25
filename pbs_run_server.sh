#!/bin/sh

#PBS -l nodes=2:ppn=8
#PBS -l walltime=96:00:00
#PBS -l pmem=16gb
cd $PBS_O_WORKDIR

module load openmpi/gnu
main_dir=${PWD}

# download data
cd data && git clone -o 20news_bydate https://github.com/bobye/20newsgroups.git 
cd 20news_bydate/20newsgroups_clean
bunzip2 20newsgroups.d2s.bz2
bunzip2 20newsgroups.d2s.vocab0.bz2


# run d2 clustering
cd $main_dir
num_of_nodes=16
data_files=data/20news_bydate/20newsgroups_clean/20newsgroups.d2s

# make sure split data only once
batch_size=`./d2 -i $data_files -n 20000 -d 300 -s 128 --type 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`


# run with different cluster numbers, use nohup to forward stdout
# results will save to different files with multiple runs
k=10 # 20 30 40 60 80
s=64 # 8 16 32
nohup mpirun -n $num_of_nodes ./d2 -i $data_files -n $batch_size -d 300 -s $s --type 7 --clusters $k > nohup.$k.out 

# save nohup.*.out for convergence analysis

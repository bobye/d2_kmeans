#!/bin/bash

echo 'Starting the test of parallel version'

#alias mpirun="mpirun -mca btl ^openib" # suppress openib

cd ..
#cd .. && make clean && make MPI=1 &> /dev/null

num_of_nodes=8
batch_size=`./d2 -i data/icip14_data/total.d2 -p 2 -n 5000 -d 3,3 -s 6,8 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`

echo "[$num_of_nodes processes] Profiling image centroids: test single phase #0 with one cluster"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 1

echo "[$num_of_nodes processes] Profiling image centroids: test single phase #1 with one cluster"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 1 --clusters 1

echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with #0 phase"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 10 --max_iter 20


echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with #1 phase"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 1 --clusters 10 --max_iter 20

echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with both phases"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --clusters 10 --max_iter 20

cd data/protein_seq
make clean & make MPI=1 

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../

./protein protein -1 -$num_of_nodes 0

echo "[$num_of_nodes processes] Profile protein centroids: test 3-gram with one cluster"
time mpirun -n $num_of_nodes ./protein protein 2 1 0

echo "[$num_of_nodes processes] Profile protein centroids: test 1,2,3-gram with one cluster"
time mpirun -n $num_of_nodes ./protein protein -1 1 0

echo "[$num_of_nodes processes] Clustering proteins: n=10742 1-gram with 2 cluster"
time mpirun -n $num_of_nodes ./protein protein 0 2 0


#num_of_nodes=4
#echo "[$num_of_nodes processes] Clustering documents: n=1000 and k=3 with both phases"
#batch_size=`./d2 -i data/20news-bydate/20newsgroups_clean/20newsgroups.d2s -n 1000 -d 300 -s 128 --type 7 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`

#time mpirun -n $num_of_nodes ./d2 -i data/20news-bydate/20newsgroups_clean/20newsgroups.d2s -n $batch_size -d 300 -s 128 --clusters 3  --type 7

#echo 'End the test of parallel version [Success]'


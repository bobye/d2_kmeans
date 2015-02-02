#!/bin/sh

echo 'Starting the test of parallel version'

alias mpirun="mpirun -mca btl ^openib" # suppress openib

cd .. && make clean && make MPI=1 &> /dev/null

echo 'Build done'

num_of_nodes=8
batch_size=`./d2 -i data/icip14_data/total.d2 -p 2 -n 5000 -d 3,3 -s 6,8 --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`

echo "[$num_of_nodes processes] Profiling image centroids: test single phase #0 with one cluster"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 1 -o centroids.d2 > /dev/null


echo "[$num_of_nodes processes] Profiling image centroids: test single phase #1 with one cluster"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 1 --clusters 1 -o centroids.d2 > /dev/null


echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with #0 phase"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 0 --clusters 10 --max_iter 20 -o centroids.d2 > /dev/null


echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with #1 phase"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --phase_only 1 --clusters 10 --max_iter 20 -o centroids.d2 > /dev/null


echo "[$num_of_nodes processes] Clustering images: n=2000 and k=10 with both phases"
time mpirun -n $num_of_nodes ./d2 -i data/icip14_data/total.d2 -p 2 -n $batch_size -d 3,3 -s 8,8 --clusters 10 --max_iter 20 -o centroids.d2 > /dev/null

cd data/protein_seq

./protein protein -1 -$num_of_nodes 0 > /dev/null

echo "[$num_of_nodes processes] Profile protein centroids: test 3-gram with one cluster'
time mpirun -n $num_of_nodes ./protein protein 2 1 0 > /dev/null

echo "[$num_of_nodes processes] Profile protein centroids: test 1,2,3-gram with one cluster"
time mpirun -n $num_of_nodes ./protein protein -1 1 0 > /dev/null

echo "[$num_of_nodes processes] Clustering proteins: n=10742 1-gram with 2 cluster"
time mpirun -n $num_of_nodes ./protein protein 0 2 0

echo 'End the test of parallel version [Success]'

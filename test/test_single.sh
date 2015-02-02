#!/bin/sh


echo 'Starting the test of serial version'

cd .. && make clean && make MPI=0 &> /dev/null

echo 'Build done'

echo 'Profiling image centroids: test single phase #0 with one cluster'
time ./d2 -i data/mountaindat.d2 -p 2 -n 1000 -d 3,3 -s 6,8 --phase_only 0 --clusters 1 -o centroids.d2 > /dev/null


echo 'Profiling image centroids: test single phase #1 with one cluster'
time ./d2 -i data/mountaindat.d2 -p 2 -n 1000 -d 3,3 -s 6,8 --phase_only 1 --clusters 1 -o centroids.d2 > /dev/null


echo 'Clustering images: n=2000 and k=10 with #0 phase'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --phase_only 0 --clusters 10 -o centroids.d2 > /dev/null


echo 'Clustering images: n=2000 and k=10 with #1 phase'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --phase_only 1 --clusters 10 -o centroids.d2 > /dev/null


echo 'Clustering images: n=2000 and k=10 with both phases'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --clusters 10 -o centroids.d2 > /dev/null

cd data/protein_seq

echo 'Profile protein centroids: test 3-gram with one cluster'
time ./protein protein 2 1 0 > /dev/null

echo 'Profile protein centroids: test 1,2,3-gram with one cluster'
time ./protein protein -1 1 0 > /dev/null

echo 'Clustering proteins: n=10742 1-gram with 2 cluster, it takes tens of seconds'
time ./protein protein 0 2 0

echo 'End the test of serial version [Success]'




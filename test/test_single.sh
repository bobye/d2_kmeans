#!/bin/bash


echo 'Starting the test of serial version'

cd ..
#cd .. && make clean && make MPI=0 &> /dev/null

echo 'Profiling image centroids: test single phase #0 with one cluster'
time ./d2 -i data/mountaindat.d2 -p 2 -n 1000 -d 3,3 -s 6,8 --phase_only 0 --clusters 1  > /dev/null


echo 'Profiling image centroids: test single phase #1 with one cluster'
time ./d2 -i data/mountaindat.d2 -p 2 -n 1000 -d 3,3 -s 6,8 --phase_only 1 --clusters 1  > /dev/null


echo 'Clustering images: n=2000 and k=10 with #0 phase'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --phase_only 0 --clusters 10  > /dev/null


echo 'Clustering images: n=2000 and k=10 with #1 phase'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --phase_only 1 --clusters 10  > /dev/null


echo 'Clustering images: n=2000 and k=10 with both phases'
time ./d2 -i data/mountaindat.d2 -p 2 -n 2000 -d 3,3 -s 6,8 --clusters 10  > /dev/null

cd data/protein_seq
make clean && make MPI=0 &> /dev/null

echo 'Profile protein centroids: test 3-gram with one cluster'
time ./protein protein 2 1 0 > /dev/null

echo 'Profile protein centroids: test 1,2,3-gram with one cluster'
time ./protein protein -1 1 0 > /dev/null

echo 'Clustering proteins: n=10742 1-gram with 2 cluster, it takes tens of seconds'
time ./protein protein 0 2 0


echo 'Profiling bag of word vectors: n=100 k=1'
time ./d2 -i data/20news-bydate/20newsgroups_clean/20newsgroups.d2s -n 100 -d 300 -s 128 --clusters 1 --type 7 > /dev/null

echo 'End the test of serial version [Success]'


#!/bin/sh


echo 'Profiling centroids: test single phase #0 with one cluster'
echo
time ./d2 -i mountaindat.txt -p 2 -n 1000 -d 3,3 -s 6,11 --phase_only 0 --clusters 1


echo 'Profiling centroids: test single phase #1 with one cluster'
echo
time ./d2 -i mountaindat.txt -p 2 -n 1000 -d 3,3 -s 6,11 --phase_only 1 --clusters 1


echo 'Clustering: n=5000 and k=10 with #0 phase'
echo
time ./d2 -i mountaindat.txt -p 2 -n 2000 -d 3,3 -s 6,11 --phase_only 0 --clusters 10

#!/bin/sh

## Parameters Set
num_of_nodes=8 # MPI processors
meta_class=50 # variability of population

dim_of_vec=64
size_of_bag=64
k=16
total_size=1280
max_iter=20

filename="${total_size}_${dim_of_vec}_${size_of_bag}_${meta_class}"


## Generate data, it takes from several seconds to minutes
cd ../matlab
octave --silent --eval "syntheticdata($dim_of_vec, $size_of_bag, $meta_class, $total_size)"
echo 'Data generated'

cd .. && make clean && make MPI=1 &> /dev/null
echo 'Build done'

## Split data into segments at the same number of processors
mkdir -p data/synthetic_data
batch_size=`./d2 -i data/synthetic_data/${total_size}_${dim_of_vec}_${size_of_bag}_${meta_class}.d2 -p 1 -n ${total_size} -d ${dim_of_vec} -s ${size_of_bag} --prepare_batches $num_of_nodes | grep batch_size | sed 's/^.*batch_size://g'`
echo "Split data with segment size $batch_size over $num_of_nodes processors"

## run the program with fixed iterations
time mpirun -n $num_of_nodes ./d2 -i data/synthetic_data/${filename}.d2 -p 1 -n ${batch_size} -d ${dim_of_vec} -s ${size_of_bag} --clusters $k --max_iter $max_iter -o data/synthetic_data/${filename}_p.d2


## The following command is the same as before 
## but no triangle inequality is used to prone 

# time mpirun -n $num_of_nodes ./d2 -i data/synthetic_data/${filename}.d2 -p 1 -n ${batch_size} -d ${dim_of_vec} -s ${size_of_bag} --clusters $k --max_iter $max_iter -o data/synthetic_data/${filename}_p.d2 --non_triangle

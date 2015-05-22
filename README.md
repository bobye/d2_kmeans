Accelerated Discrete Distributions Clustering under Wasserstein Distance
=============

AD2-clustering is an accelerated clustering algorithm for **D**iscrete **D**istributions
under the *exact* [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It can scale to large-scale D2 data with parallel computing. Please see our paper for technical details. 

[Jianbo Ye](http://www.personal.psu.edu/jxy198), [Panruo Wu](http://www.cs.ucr.edu/~pwu011/), James Z. Wang and Jia Li, "Fast Discrete Distribution Clustering under Wasserstein Distance" in submission to NIPS 2015.

This project is the __first__ efficient implementation public available for D2-clustering,
and it is still in the very early stage of development. Hence please don't expect any
guanrantee for backward compatibility at this point. 

## Discrete Distributions
There are three data types of discrete distribution that have been covered
in this project:
 - [default] discrete distribution over vector space endowed with Euclidean distance
 - Same type as the default one, but over embeddings space with a finite vocabulary size
   (e.g. document represented by tf-idf over word vectors, sparsified histograms)
 - normalized (dense) histograms with bin-to-bin distance
 - [**experimental] d2 over [n-gram](http://en.wikipedia.org/wiki/N-gram) provided
   with item-to-item similarity/distance (sparse histograms are represented as 1-gram)

An object/instance can be represented as the joint of multiple discrete
distributions (of any aforementioned type), called phases. For example, an image can be
represented in two phases: color distribution and texture distribution; a protein
sequence can be represented as distributions in three phases, aka,
1-gram, 2-gram, 3-gram of amino acid. A document can be represented by 
a bag of weighted word vectors.
The co-clustering is then performed jointly over multiple phases.

## How to compile and run tests
Though one can build a serial version to solve small-to-moderate scale problems with speed,
the strength of AD2-clustering is its scalability to large dataset with parallelization at 
high scaling efficiency (over 80% on hundreds of cores). 

Build dependencies:
 - MPI
 - CBLAS
 - [Mosek](https://mosek.com) 7.0+: free individual academic license available.


Create `make.inc` file to resolve the dependencies, see `make.inc.Linux` for an example.
One can choose to build the program under Linux or Mac. 

Compile
```
 $ make MPI=0 # build sequential version, or
 $ make MPI=1 # build MPI version (default)
```

Run unit tests (it takes several minutes):
```
 $ make test 
```

By default, 64bit floating numbers are used. Optionally, you can switch to use 32bit to save
memory and get some marginal speedups by changing the line in `include/common.h`
```c
#define _D2_DOUBLE // to _D2_SINGLE
```

## Usage

### Arguments
Input options
 - `--ifile <input_filename>, -i <input_filename>` : the main name of input file (required).
 - `-n <integer>` : number of instances read in per processor (required)
 - `--phase <integer>, -p <integer>` : the number of phases per instance in `<input_filename>` (default: 1).
 - `--phase_only <integer>, -t <integer>` : the phase to cluster upon (default: use all phases)
 - `--strides <integer array>, -s <integer array>` : the numbers of support points of computed centroids in each phase (required), integer array with comma delimiter and no spaces. 
 - `-d <integer array>` : the dimensions in each phase (required), integer array with comma delimiter and no spaces.
 - `--types <integer>, -E <integer>` : the type of D2 data (default: 0, see `include/d2_param.h` for details).
 
Algorithm options
 - `--clusters <integer>, -c <integer>` : number of clusters intended (default: 3, mostly required).
 - `--max_iters <integer>, -m <integer>` : the maximal number of iterations (default: 100).
 - `--non_triangle, -T` : disable the triangle inequality based acceleration (default: enabled).
 - `--eval <centroids_filename>, -e <centroids_filename>` : no clustering, but assigning instances to the pre-computed centroids (default: disabled).
 
Parallel computing options
 - `--prepare_batches <integer>, -P <integer>` : the number of batches that is equal to the number of processors in data pre-processing stage. It reads in `<input_filename> = data.d2` and generates files in say `data.d2.0, data.d2.1, data.d2.2, data.d2.3` containing randomly splitted parts of `data.d2`.
 
### Modes
1. The default mode is the clustering algorithm, which outputs results in two main files. For example, taking in `data.d2`, program outputs the centroids computed (`data.d2_123456_c.d2` which is again in D2 format) and the memberships of each instances (`data.d2_123456.label` in sequential run or `data.d2_123456.label_o` in parallel run). 
2. To preprocess a data file into multiple batches and later feed them into a parallel computing environment, one has to call `--prepare_batches`.
3. Given pre-computed centroids from a training set, one can assign cluster memberships to another testing set using `--eval`. 

### Data formats
See [data](data) for options to prepare D2 data. 

## Examples
 - [preprocessed "20newsgroups" as bags of word-vectors](pbs_run_server.sh) (related [codes](https://github.com/bobye/20newsgroups) on preparing the data).
 - [preprocessed n-gram of protein sequences](data/protein_seq): experimental
 - [USPS handwritten digit dataset](pbs_run_server_usps.sh)

See [test](test) for more examples.




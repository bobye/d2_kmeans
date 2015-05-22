Accelerated Discrete Distributions Clustering under Wasserstein Distance
=============

AD2-clustering is an accelerated clustering algorithm for **D**iscrete **D**istributions (D2)
under the *exact* [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It can scale to large-scale D2 data with parallel computing. Please see our paper for technical details: 

[Jianbo Ye](http://www.personal.psu.edu/jxy198), Panruo Wu, James Z. Wang and [Jia Li](http://sites.stat.psu.edu/~jiali/), "Fast Discrete Distribution Clustering under Wasserstein Distance" in submission to NIPS 2015.

This project is the __first__ efficient implementation public available for D2-clustering,
and it is still in the very early stage of development.

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

See [here](data) for options to prepare D2 data and 
see [here](src) for detailed instructions on possible arguments and modes of program.

## Examples
 - [preprocessed "20newsgroups" as bags of word-vectors](experiment/pbs_run_server_20news.sh) (related [codes](https://github.com/bobye/20newsgroups) on preparing the data).
 - [preprocessed n-gram of protein sequences](data/protein_seq): experimental
 - [USPS handwritten digit dataset](experiment/pbs_run_server_usps.sh)

See [here](test) for more examples.




Fast Discrete Distribution Clustering Using Wasserstein Barycenter with Sparse Support
=============

AD2-clustering is a fast clustering algorithm for **D**iscrete **D**istributions (D2)
under the *exact* [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It can scale to large-scale D2 data with parallel computing. 
This project is the __first__ efficient implementation publicly available for D2-clustering,
and it is still in the very early stage of development. Contributions are welcomed. 

Please see our paper for technical details: 

[Jianbo Ye](http://www.personal.psu.edu/jxy198), [Panruo Wu](http://www.cs.ucr.edu/~pwu011/), James Z. Wang and Jia Li,
Fast Discrete Distribution Clustering Using Wasserstein Barycenter with Sparse Support,
IEEE Transactions on Signal Processing, 2017 ([arXiv:1510.00012](http://arxiv.org/abs/1510.00012) [stat.CO], September 2015)


## Discrete Distributions
There are four data types of discrete distribution that have been covered
in this project (See [data format](data) for their specifications):
 - [default] discrete distribution over vector space endowed with Euclidean distance
 - Same type as the default one, but over embeddings space with a finite vocabulary size
   (e.g. document represented by tf-idf over word vectors, sparsified histograms)
 - normalized dense histograms with bin-to-bin distance
 - normalized sparse histograms with bin-to-bin distance
 - [**experimental] d2 over [n-gram](http://en.wikipedia.org/wiki/N-gram) provided
   with item-to-item similarity/distance 
   (sparse histograms as the centroids are represented as 1-gram)

An object/instance can be represented as the joint of multiple discrete
distributions (in combinations of any aforementioned types), called phases. 
For example, an image can be represented in two phases: 
color distribution and texture distribution; a protein
sequence can be represented as distributions in three phases, aka,
1-gram, 2-gram, 3-gram of amino acid. A document can be represented by 
a bag of weighted word vectors.
The co-clustering is then performed jointly over multiple phases.

## How to compile and run tests
Though one can build a serial version to solve small-to-moderate scale problems with speed,
the strength of AD2-clustering is its scalability to large dataset with parallelization at 
high scaling efficiency (over 80% on hundreds of cores). You can find [quickstart tips on 
the Google cloud compute engine](https://github.com/bobye/d2_kmeans/wiki) for references.

Build dependencies:
 - MPI
 - CBLAS
 - [Mosek](https://mosek.com) 7.x: free individual academic license available.


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
memory in computation and get some marginal speedups by edit the line in `make.inc`
```makefile
D2_DEFINES=-D _D2_DOUBLE # change to _D2_SINGLE if the size of RAM is limited
```

## Usage

See [here](data) for options to prepare D2 data and 
see [here](src/app) for detailed instructions on possible arguments and modes of program.

## Examples
 - [preprocessed "20newsgroups" as bags of word-vectors](experiment/pbs_run_server_20news.sh) (related [codes](https://github.com/bobye/20newsgroups) on preparing the data).
 - [preprocessed n-gram of protein sequences](data/protein_seq): experimental
 - [USPS handwritten digit dataset](experiment/pbs_run_server_usps.sh)

See [here](test) for more examples.


## Misc
 - Related projects: [d2suite](https://github.com/colourbrain/d2suite)



d2_kmeans: accelerated and parallel clustering of discrete distributions
=============

D2-Kmeans is a clustering algorithm for __D__iscrete __D__istributions,
including normalized histogram and correlated n-gram as special cases,
under the [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It follows the Elkan's algorithm [1] of k-means to save unnecessary
computation of instance re-labeling, and uses a modifed version of
Bregman ADMM [2] to accelerate and scale-up the computation of cluster centroids.
Please see [3] for technical details. 

## Discrete Distributions
There are three data types of discrete distribution that have been covered
in this project:
 - (*default*) discrete distribution over vector space endowed with Euclidean distance
 - Same type as the default one, but with d2 over word embeddings
   (e.g. document represented by tf-idf over word vectors)
 - normalized (dense) histograms with bin-to-bin distance
 - d2 over [n-gram](http://en.wikipedia.org/wiki/N-gram) provided
   with item-to-item similarity/distance (sparse histograms are represented as 1-gram)

An object/instance can be represented as the joint of multiple discrete
distributions (of any aforementioned type), called phases. For example, an image can be
represented in two phases: color distribution and texture distribution; a protein
sequence can be represented as distributions in three phases, aka,
1-gram,2-gram,3-gram of amino acid.
The co-clustering is then performed jointly over multiple phases.

## How to compile and run tests

Build dependencies:
 - MPI
 - CBLAS: OpenBLAS or MKL or Accelerate(Mac)
 - [Mosek](https://mosek.com) 7.0+
 
```
 $ make MPI=0 # build sequential version, or
 $ make MPI=1 # build MPI version
```

Run unit tests:
```
 $ make test # it takes several minutes. 
```

## Guides and Tutorials
 - [data format supported](data)
 - [clustering "20newsgroups" as bags of word-vectors](https://github.com/bobye/20newsgroups/wiki)

## Reference
1. Elkan, Charles. "Using the triangle inequality to accelerate k-means." ICML. Vol. 3. 2003.
2. Wang, Huahua et al. "Bregman Alternating Direction Method of Multipliers." NIPS. 2014.
3. Ye, Jianbo et al. "AD2-Clustering: Accelerated Clustering of Discrete Distributions with Bregman ADMM" in submission.




d2_kmeans: accelerated and parallel clustering of discrete distributions
=============

D2-Kmeans is a clustering algorithm for discrete distributions,
including normalized histogram as a special case,
under the [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It follows the Elkan's algorithm[1] of k-means to save unnecessary
computation of instance re-labeling, and uses Bregman ADMM[2] to accelerate
and scale-up the computation of cluster centroids.
Please see [3] for technical details. 

## Discrete Distributions
There are three data types of discrete distribution that have been covered
in this project:
 - discrete distribution over vector space endowed with Euclidean distance
 - normalized (dense) histograms with bin-to-bin distance
 - discrete distribution over [n-gram](http://en.wikipedia.org/wiki/N-gram)
   with item-to-item similarity/distance (sparse histograms are represented
   as 1-gram)

An object/instance can be represented as the joint of multiple discrete
distributions (of any aforementioned type). For example, an image can be
represented as color distribution and texture distribution; a protein
sequence can be represented as distributions in three phases, aka,
1-gram,2-gram,3-gram of amino acid.
The clustering is then performed jointly over differernt phases.

## How to compile and run tests

dependencies:
 - MPI
 - CBLAS: OpenBLAS or MKL or Accelerate(Mac)
 - [Mosek](https://mosek.com) 7.0+
 
```
 $ make MPI=0 # build sequential version, or
 $ make MPI=1 # build MPI version
```

tests (TBA):
```
 $ make test # it takes several minutes. 
```

## Guides and Tutorials
 - a brief user guide
 - [clustering documents as bags of word-vectors](https://github.com/bobye/d2_kmeans/wiki/Document-Clustering)

## Port to Interfaces
 - Matlab [to appear]

## Reference
1. Elkan, Charles. "Using the triangle inequality to accelerate k-means." ICML. Vol. 3. 2003.
2. Wang, Huahua et al. "Bregman Alternating Direction Method of Multipliers." NIPS. 2014.
3. Ye, Jianbo et al. "AD2-Clustering: Accelerated Clustering of Discrete Distributions with Bregman ADMM" in submission.

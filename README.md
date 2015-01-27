d2_kmeans: accelerated and parallel clustering of discrete distributions
=============

D2-Kmeans is a clustering algorithm for discrete distributions,
including normalized histogram as a special case,
under the [Wasserstein metric](http://en.wikipedia.org/wiki/Wasserstein_metric).
It follows the Elkan's algorithm[1] of k-means to save unnecessary
computation of instance re-labeling, and uses Bregman ADMM[2] to accelerate
and scale-up the computation of cluster centroids.
Please see [3] for technical details. 

## How to compile

### Linux + Intel Compiler
Parallel version is based on Linux and Intel Compiler with dependencies

 - MKL
 - openMP support
 - Mosek 7.0+

### Mac OSX 10.10 + Clang
Besides, a serial version is available for Mac OSX whose dependencies include

 - Accelerate Framework
 - no OpenMP support
 - Mosek 7.0+

## Interfaces
 - Matlab [to appear]


## Reference
1. Elkan, Charles. "Using the triangle inequality to accelerate k-means." ICML. Vol. 3. 2003.
2. Wang, Huahua et al. "Bregman Alternating Direction Method of Multipliers." NIPS. 2014.
3. Ye, Jianbo et al. "D2-Kmeans: Accelerated and Parallel Clustering of Discrete Distributions" in submission.

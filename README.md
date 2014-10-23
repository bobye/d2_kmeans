d2_clustering
=============

scalable and parallel D2 (discrete distribution) clustering

- Algorithm Prototypes
```
  matlab/
    d2_clustering.m --- old scripts
    d2clusters.m --- old clustering framework
    profile_centroid.m --- profiling the convergence of centroid updates
    profile_kantorovich.m --- profiling LP solution of transportation problem
    centroid_sph*.m --- computing centroid of a single phase
```
- C/C++ sources (algorithms should be implemented in C99 std, but tests or others can be in C++)
```
  include/
  src/
  Makefile -- OSX Accelerate (for BLAS and LAPACK)
  Makefile.MKL --- ICC/ICPC and MKL (for BLAS and LAPACK)
```
- Dataset
We implement multi-phases D2, which means each object can be represented
by multiple D2. For example, an image can be represented by a D2 in color
space and a D2 in texture space. To read a 2-phase 1000 entries with the 
first phase in 3 dimension and 6 avg bins, and the second phase  in 3 dimension
and 11 bins. You may type
```bash
time ./d2 mountaindat.txt 2 1000 3,3 6,11
```

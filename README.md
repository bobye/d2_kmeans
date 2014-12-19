d2_clustering
=============

scalable and parallel D2 (discrete distribution) clustering

## Dataset and Format
We implement multi-phases D2, which means each object can be represented
by multiple D2. For example, an image can be represented by a D2 in color
space and a D2 in texture space. The format of data is as follows:
```emacs-lisp
;; start of first object
d1 ;; dimension of the first phase
n1 ;; number of bins in the first phase
w{1} w{2} ... w{n1} ;; weights of bins
x{1,1} x{1,2} ... x{1,d1}
x{2,1} x{2,2} ... x{2,d1}
...
x{n1,1} x{n1,2} ... x{n1,d1}
d2 ;; dimension of the second phase
n2 ;; number of bins in the second phase
w{1} w{2} ... w{n2} 
x{1,1} x{1,2} ... x{1,d2}
x{2,1} x{2,2} ... x{2,d2}
...
x{n2,1} x{n2,2} ... x{n2,d2}
;; end of first object
;; start of the second object
...
```
It is required that w{i} are strictly larger than zero.

To read a 2-phase 1000 entries with the first phase in 3 dimension
and 6 average number of bins, and the second phase  in 3 dimension
and 11 bins. You may type
```bash
time ./d2 mountaindat.txt 2 1000 3,3 6,11
```


## Algorithm Prototypes
It contains prototype version in Matlab.
```
  matlab/
    d2_clustering.m --- old scripts
    d2clusters.m --- old clustering framework
    profile_centroid.m --- profiling the convergence of centroid updates
    profile_kantorovich.m --- profiling LP solution of transportation problem
    centroid_sph*.m --- computing centroid of a single phase
```

## C/C++ sources
Core algorithms, aka main source, should be implemented in C99 std, but tests or others can be in C++ or other languages.
```
  include/
  src/
  Makefile -- OSX Accelerate (for BLAS and LAPACK)
  Makefile.MKL --- ICC/ICPC and MKL (for BLAS and LAPACK)
```

d2_kmeans
=============

D2-Kmeans: accelerated and parallel clustering of discrete distributions

## How to compile

### Linux + Intel Compiler

dependencies

 - MKL
 - openMP support
 - Mosek 7.0+

### Mac OS 10.10 + Clang

dependencies

 - Accelerate Framework
 - no OpenMP support
 - Mosek 7.0+


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
It is required that `w{i}` are strictly larger than zero.

To read 1000 2-phase entries with the first phase in 3 dimension
and 6 average number of bins, and the second phase  in 3 dimension
and 11 bins. You may type
```bash
time ./d2 -i mountaindat.txt -p 2 -n 1000 -d 3,3 -s 6,11
```

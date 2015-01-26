## Dataset and Format
We provides two kinds of input formats

1. We implement multi-phases D2, which means each object can be represented
   by multiple D2. For example, an image can be represented by a D2 in color
   space and a D2 in texture space. The format of .d2 data is as follows:
```emacs-lisp
;; filen ext: .d2
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
2. In some cases, it would be useful to work with histograms (where the sum
   of bins is equal to one). The format of data is as follows:

```emacs-lisp
;; file ext: .hist.d
n ;; number of bins
d{1,1} d{2,1} ... d{n,1} ;; cost of transportation between different bins
d{1,2} d{2,2} ... d{n,2} ;; d has to be symmetric, and d{i,i} = 0
...
d{1,n} d{2,n} ... d{n,n}
```
   Generally, it is not necessary that the distance d^{1/p} is a true metric:
   But for using triangle inequality to accelerate the undergoing computation,
   you can enforce to modify the distance such that they are qualified under
   a true metric. 

#include "blas_like.h"


// a(:,*) = a(:,*) .+ b
void _dgcmv(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pa += *pb;
}

// a(*,:) = a(*,:) .+ b
void _dgrmv(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa += *pb;
}

// a = diag(b) * a
void _dgcms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pa *= *pb;
}

// a = a * diag(b) 
void _dgrms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa += *pb;
}

// b(*) = sum(a(:,*))
void _dcsum(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a,pa=b; i<n; ++i, ++pb) {
    *pb = 0;
    for (j=0; j<m; ++j, ++pa)
      *pb += *pa;
  }
}


// b(*) = sum(a(*,:))
void _drsum(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (j=0,pb=b; j<m; ++j, ++pb) 
    *pb = 0;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pb += *pa;  
}

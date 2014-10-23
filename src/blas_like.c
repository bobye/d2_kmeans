#include <math.h>

#define __USE_C99_MATH
#include <stdbool.h>
#include "blas_like.h"
#include <stdio.h>

#define SCALAR double

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_FREE(x)         free(x)

#elif defined __USE_MKL__
#include <mkl.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) mkl_malloc( (x) *sizeof(SCALAR), 16) 
#define _D2_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _D2_FREE(x)         mkl_free(x)

#elif defined __GNUC__
#include <cblas.h>
#include <lapacke.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_FREE(x)         free(x)

#endif


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
  double *pa =a, *pb =b;
  for (i=0; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa += *pb;
}

// a = diag(b) * a
void _dgcms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa = a, *pb;
  for (i=0; i<n; ++i)
    for (j=0, pb=b; j<m; ++j, ++pa, ++pb)
      *pa *= *pb;
}

// a = a * diag(b) 
void _dgrms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa = a, *pb = b;
  for (i=0; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa *= *pb;
}

// a = diag(1./b) * a
void _dicms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a; i<n; ++i)
    for (j=0,pb=b; j<m; ++j, ++pa, ++pb)
      *pa /= *pb;
}

// a = a * diag(1./b) 
void _dirms(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i,++pb)
    for (j=0; j<m; ++j, ++pa)
      *pa /= *pb;
}

// b(*) = sum(a(:,*))
void _dcsum(int m, int n, double *a, double *b) {
  int i,j;
  double *pa, *pb;
  for (i=0,pa=a,pb=b; i<n; ++i, ++pb) {
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

// normalize by column
void _dcnorm(int m, int n, double *a, double *sa) {
  int i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _D2_MALLOC_SCALAR(n);
  }    
  _dcsum(m, n, a, sa);
  for (i=0, pa=sa; i<n; ++i, ++pa) {
    for (j=0; j<m; ++j, ++a) (*a) /= *pa;
  }
  if (!isAllocated) _D2_FREE(sa);
}

// normalize by row
void _drnorm(int m, int n, double *a, double *sa) {
  int i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = _D2_MALLOC_SCALAR(m);
  }    
  _drsum(m, n, a, sa);
  for (i=0; i<n; ++i) {
    pa = sa;
    for (j=0; j<m; ++j, ++a, ++pa) (*a) /= *pa;
  }
  if (!isAllocated) _D2_FREE(sa);
}

// c = a.*b
void _dvmul(int n, double *a, double *b, double *c) {
  int i;
  for (i=0; i<n; ++i, ++c, ++a, ++b)
    *c = (*a) * (*b);
}

void _dpdist2(int d, int n, int m, double * A, double * B, double *C) {
  double *AA, *BB, *sA, *sB;
  //assert(d>0 && n>0 && m>0);

  AA = _D2_MALLOC_SCALAR(d*n);
  BB = _D2_MALLOC_SCALAR(d*m);
  sA = _D2_MALLOC_SCALAR(n);
  sB = _D2_MALLOC_SCALAR(m);

  _dvmul(d*n, A, A, AA);
  _dvmul(d*m, B, B, BB);

  _dcsum(d, n, AA, sA);
  _dcsum(d, m, BB, sB);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, m, d, -2, 
	      A, d, B, d, 0, C, n);
  
  _dgcmv(n, m, C, sA);
  _dgrmv(n, m, C, sB);
  
  _D2_FREE(AA);
  _D2_FREE(BB);
  _D2_FREE(sA);
  _D2_FREE(sB);
}

// inplace a -> exp(a)
void _dexp(int n, double *a) {
  int i;
  for (i=0; i<n; ++i, ++a) *a = exp(*a);
}

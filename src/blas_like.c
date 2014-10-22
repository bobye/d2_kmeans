#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

#elif defined __INTEL_COMPILER

#include <mkl.h>

#elif defined __GNUC__
#include <math.h>
#include <cblas.h>

#endif

#define __USE_C99_MATH
#include <stdbool.h>

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
    sa = (double *) malloc(n *sizeof(double) );
  }    
  _dcsum(m, n, a, sa);
  for (i=0, pa=sa; i<n; ++i, ++pa) {
    for (j=0; j<m; ++j, ++a) (*a) /= *pa;
  }
  if (!isAllocated) free(sa);
}

// normalize by row
void _drnorm(int m, int n, double *a, double *sa) {
  int i, j;
  double *pa;
  bool isAllocated = true;
  if (!sa) {
    isAllocated = false;
    sa = (double *) malloc(m *sizeof(double) );
  }    
  _drsum(m, n, a, sa);
  for (i=0; i<n; ++i) {
    pa = sa;
    for (j=0; j<m; ++j, ++a, ++pa) (*a) /= *pa;
  }
  if (!isAllocated) free(sa);
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

  AA = (double*) malloc(d*n*sizeof(double));
  BB = (double*) malloc(d*m*sizeof(double));
  sA = (double*) malloc(n*sizeof(double));
  sB = (double*) malloc(m*sizeof(double));

  _dvmul(d*n, A, A, AA);
  _dvmul(d*m, B, B, BB);

  _dcsum(d, n, AA, sA);
  _dcsum(d, m, BB, sB);

  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, n, m, d, -2, 
	      A, d, B, d, 0, C, n);
  
  _dgcmv(n, m, C, sA);
  _dgrmv(n, m, C, sB);
  
  free(AA);
  free(BB);
  free(sA);
  free(sB);
}

// inplace a -> exp(a)
void _dexp(int n, double *a) {
  int i;
  for (i=0; i<n; ++i, ++a) *a = exp(*a);
}

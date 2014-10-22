#ifndef _BLAS_LIKE_H_
#define _BLAS_LIKE_H_


#ifdef __cplusplus
extern "C" {
#endif

  #include <stdlib.h>

  // element-wise op
  void _dvmul(int n, double *a, double *b, double *c);// c = a .* b
  void _dexp(int n, double *a);//inplace a -> exp(a);

  // column-wise op
  void _dgcmv(int m, int n, double *a, double *b); // a(:,*) = a(:,*) .+ b
  void _dgcms(int m, int n, double *a, double *b); // a = diag(b) * a
  void _dicms(int m, int n, double *a, double *b); // a = diag(1./b) * a
  void _dcsum(int m, int n, double *a, double *b); // b(*) = sum(a(:,*))
  void _dcnorm(int m, int n, double *a, double *sa); // inplace a(:,*) -> a(:,*) / sum(a(:,*))

  // row-wise op
  void _dgrmv(int m, int n, double *a, double *b); // a(*,:) = a(*,:) .+ b
  void _dgrms(int m, int n, double *a, double *b); // a = a * diag(b) 
  void _dirms(int m, int n, double *a, double *b); // a = a * diag(1./b) 
  void _drsum(int m, int n, double *a, double *b); // b(*) = sum(a(*,:))
  void _drnorm(int m, int n, double *a, double *sa); // inplace a(*,:) = a(*,:) / sum(a(*,:))


  /* compute squared Euclidean distance matrix
   * A: d x n 
   * B: d x m
   * C: n x m 
   * d: dimension of data entry
   * n, m: number of data entry
   */
  void _dpdist2(int d, int n, int m, double * A, double * B, double *C);
  

#ifdef __cplusplus
}
#endif


#endif /* _BLAS_LIKE_H_ */

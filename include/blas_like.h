#ifndef _BLAS_LIKE_H_
#define _BLAS_LIKE_H_


#ifdef __cplusplus
extern "C" {
#endif

  #include <stdlib.h>

  // assertation
  void _dgzero(size_t n, double *a); //assert (a>0)

  // element-wise op
  void _dvmul(size_t n, double *a, double *b, double *c);// c = a .* b
  void _dexp(size_t n, double *a);//inplace a -> exp(a);

  // column-wise op
  void _dgcmv(size_t m, size_t n, double *a, double *b); // a(:,*) = a(:,*) .+ b
  void _dgcms(size_t m, size_t n, double *a, double *b); // a = diag(b) * a
  void _dicms(size_t m, size_t n, double *a, double *b); // a = diag(1./b) * a
  void _dcsum(size_t m, size_t n, double *a, double *b); // b(*) = sum(a(:,*))
  void _dcsum2(size_t m, size_t n, double *a, double *b); // b(*) += sum(a(:,*))
  void _dcnorm(size_t m, size_t n, double *a, double *sa); // replace a(:,*) -> a(:,*) / sum(a(:,*))
  void _dccenter(size_t m, size_t n, double *a, double *sa); // replace a(:,*) -> a(:,*) - mean(a(:,*))
  // row-wise op
  void _dgrmv(size_t m, size_t n, double *a, double *b); // a(*,:) = a(*,:) .+ b
  void _dgrms(size_t m, size_t n, double *a, double *b); // a = a * diag(b) 
  void _dirms(size_t m, size_t n, double *a, double *b); // a = a * diag(1./b) 
  void _drsum(size_t m, size_t n, double *a, double *b); // b(*) = sum(a(*,:))
  void _drsum2(size_t m, size_t n, double *a, double *b); // b(*) += sum(a(*,:))
  void _drnorm(size_t m, size_t n, double *a, double *sa); // inplace a(*,:) = a(*,:) / sum(a(*,:))
  void _drcenter(size_t m, size_t n, double *a, double *sa); // replace a(*,:) -> a(*,:) - mean(a(*,:))


  /* compute squared Euclidean distance matrix
   * A: d x n 
   * B: d x m
   * C: n x m 
   * d: dimension of data entry
   * n, m: number of data entry
   */
  void _dpdist2(int d, size_t n, size_t m, double * A, double * B, double *C);
  void _dpdist_symbolic(int d, size_t n, size_t m, int * A, int * B, double *C, 
			const int vocab_size, const double* dist_mat);
  

#ifdef __cplusplus
}
#endif


#endif /* _BLAS_LIKE_H_ */

#ifndef _BLAS_LIKE_H_
#define _BLAS_LIKE_H_


#ifdef __cplusplus
extern "C" {
#endif

  // column-wise op
  void _dgcmv(int m, int n, double *a, double *b); // a(:,*) = a(:,*) .+ b
  void _dgcms(int m, int n, double *a, double *b); // a = diag(b) * a
  void _dcsum(int m, int n, double *a, double *b); // b(*) = sum(a(:,*))

  // row-wise op
  void _dgrmv(int m, int n, double *a, double *b); // a(*,:) = a(*,:) .+ b
  void _dgrms(int m, int n, double *a, double *b); // a = a * diag(b) 
  void _drsum(int m, int n, double *a, double *b); // b(*) = sum(a(*,:))
  

#ifdef __cplusplus
}
#endif


#endif /* _BLAS_LIKE_H_ */

#ifndef _STAT_H_
#define _STAT_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "d2_clustering.h"


#include "blas_like.h"
#include "blas_util.h"

  /*
  void d2_mean(sph * data, int * label, long num_of_entries, int num_of_labels, 
	       __OUT__ SCALAR * means, __OUT__ SCALAR * covs);
  void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int d, int n, __OUT__ SCALAR * sample);
  */
  void shuffle(size_t * array, size_t n);

  void sp_alloc(int rows, int cols, sparse_matrix *spmat, int nnz);
  void sparse(double *mat, int rows, int cols, sparse_matrix *spmat, int nnz);
  void multdense(sparse_matrix *spmat, int m, int n, int k, double *mat1, double *mat2); // mat2 (m x k) = spmat (m x n) * mat1 (n x k)

#ifdef __cplusplus
}
#endif


#endif /* _STAT_H_ */

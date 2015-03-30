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

  typedef struct {
    int m, n, nnz; /* number of rows, number of colums, # of triplet entries */
    int *p; /* column pointers (size n+1) */
    int *i; /* row indices */
    SCALAR *x; /* numerical values */
  } sparse_matrix; // matrix in compressed-column 

#ifdef __cplusplus
}
#endif


#endif /* _STAT_H_ */

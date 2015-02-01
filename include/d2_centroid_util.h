#ifndef _D2_CENTROID_UTIL_H_
#define _D2_CENTROID_UTIL_H_


#include <float.h>

/**
 * inline functions
 */
inline void accumulate_symbolic(int d, int m, int n, const int *supp /* d x n*/, const double *xx /* m x n*/, double *z /* vocab_size x d x m  */, int vocab_size) {
  int i,j,k; double val;
  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j) 
      if ((val = xx[i*m + j]) > 1E-10) {
	for (k=0; k<d; ++k)
	  z[vocab_size*(d*j + k) + supp[i*d + k]] += val;      
      }
}

inline void minimize_symbolic(int d, int m, int *supp, const double *z, const int vocab_size, const double *dist_mat, double *z_buffer) {
  int i,j, min_idx;
  double min ;
  _D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasNoTrans, 
		       vocab_size, d * m, vocab_size, 1, dist_mat, vocab_size, z, vocab_size, 0, 
		       z_buffer, vocab_size);
  
  for (i=0; i<d*m; ++i, z_buffer += vocab_size) {
    min = DBL_MAX; 
    for (j=0; j<vocab_size; ++j) {
      if (z_buffer[j] < min) {min = z_buffer[j]; min_idx = j;}
    }
    supp[i] = min_idx;
  }
}

void calculate_distmat(sph *data_ph, int* label, size_t size, sph *c, SCALAR* C);


void broadcast_centroids(mph *centroids, int i);

#endif /* _D2_CENTROID_UTIL_H_ */

#ifndef _D2_CENTROID_UTIL_H_
#define _D2_CENTROID_UTIL_H_

#ifdef __USE_MPI__
#include <mpi.h>
#endif

#include "d2/param.h"
#include "utils/blas_util.h"
#include <float.h>

/**
 * inline functions
 */
inline void accumulate_symbolic(int d, int m, int n, const int *supp /* d x n*/, const SCALAR *xx /* m x n*/, SCALAR *z /* vocab_size x d x m  */, int vocab_size) {
  int i,j,k; double val;
  for (i=0; i<n; ++i)
    for (j=0; j<m; ++j) 
      if ((val = xx[i*m + j]) > 1E-10) {
	for (k=0; k<d; ++k)
	  z[vocab_size*(d*j + k) + supp[i*d + k]] += val;      
      }
}

inline void minimize_symbolic(int d, int m, int *supp, const SCALAR *z, const int vocab_size, const SCALAR *dist_mat, SCALAR *z_buffer) {
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



inline void broadcast_centroids(mph *centroids, int i) {
#ifdef __USE_MPI__
  /* initialize centroids from one node, and broadcast to other nodes */
  switch (centroids->ph[i].metric_type) {
  case D2_EUCLIDEAN_L2:
    MPI_Allreduce(MPI_IN_PLACE, 
		  centroids->ph[i].p_supp, 
		  centroids->ph[i].col * centroids->ph[i].dim, 
		  MPI_SCALAR, 
		  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, 
		  centroids->ph[i].p_w, 
		  centroids->ph[i].col, 
		  MPI_SCALAR,
		  MPI_SUM, MPI_COMM_WORLD);
    break;
  case D2_HISTOGRAM:
    MPI_Allreduce(MPI_IN_PLACE,
		  centroids->ph[i].p_w, 
		  centroids->ph[i].col, 
		  MPI_SCALAR,
		  MPI_SUM, MPI_COMM_WORLD);
    break;
  case D2_N_GRAM:
    MPI_Allreduce(MPI_IN_PLACE,
		  centroids->ph[i].p_supp_sym, 
		  centroids->ph[i].col * centroids->ph[i].dim, 
		  MPI_INT,
		  MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,
		  centroids->ph[i].p_w, 
		  centroids->ph[i].col, 
		  MPI_SCALAR,
		  MPI_SUM, MPI_COMM_WORLD);
    break;
  }
#endif
}


void merge         (const int dim, 
		    const SCALAR * m_supp, const SCALAR * m_w, const int m, 
		    SCALAR * c_supp, SCALAR * c_w, const int n);

#endif /* _D2_CENTROID_UTIL_H_ */

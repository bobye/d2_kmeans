//#define _VERBOSE_OUTPUT

#include <math.h>
#include "d2_math.h"
#include <stdio.h>
#include <time.h>       /* time */
#include <assert.h>


void sp_alloc(int rows, int cols, sparse_matrix *spmat, int nnz) {
  spmat->m = rows;
  spmat->n = cols;
  spmat->nnz = nnz;
  spmat->p = (int*) malloc((cols+1)*sizeof(int));
  spmat->i = (int*) malloc(nnz*sizeof(int));
  spmat->x = (double*) malloc(nnz*sizeof(double));  
}
void sparse(double *mat, int rows, int cols, sparse_matrix *spmat, int nnz) { 
  /* assume variables in mat are positive */
  int i,j,count;
  count = 0;
  for (i=0; i<cols; ++i) {
    spmat->p[i] = count;
    for (j=0; j<rows; ++j) 
      if (mat[i*rows + j] > 1E-9) {
	spmat->i[count] = j;
	spmat->x[count] = mat[i*rows + j];
	count ++;
	//printf("%d %d %f\n", i, j, mat[i*rows + j]);
      }
  }
  //printf("\n");
  spmat->p[cols]=count;
  assert(count <= nnz);
}

// mat2 (m x k) = spmat.t (m x n) * mat1 (n x k)  
void multdense(sparse_matrix *spmat, int m, int n, int k, double *mat1, double *mat2) {
  int i, j, l;
  for (i=0; i<m*k; ++i) mat2[i] = 0.f;
  for (i=0; i<m; ++i) {
    for (j=spmat->p[i]; j<spmat->p[i+1]; ++j) {
      int ii = spmat->i[j];
      double xx = spmat->x[j];
      double *mat1xx = mat1 + ii;
      double *mat2xx = mat2 + i;
      for (l=0; l<k; ++l, mat1xx +=n, mat2xx += m) 
	{ *mat2xx += xx * *mat1xx; }
    }
  }
}


double randn () {
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
 
  if (call == 1) {
    call = !call;
    return X2;
  }

  do {
    U1 = -1 + ((double) rand() / RAND_MAX) * 2;
    U2 = -1 + ((double) rand() / RAND_MAX) * 2;
    W = pow(U1, 2) + pow(U2, 2);
  }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return X1;
}


/* Arrange the N elements of ARRAY in random order (no seeds).
   Only effective if N is much smaller than RAND_MAX;
   if this may not be the case, use a better random
   number generator. */
void shuffle(size_t *array, size_t n)
{
    srand (time(NULL) * world_rank);
    if (n > 1) 
    {
      int k;
      for (k=0; k<10; ++k) {
        size_t i;
        for (i = 0; i < n - 1; i++) 
        {
          size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
          size_t t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
      }
    }
}
/*
void d2_mean(sph * data, int * label, long num_of_entries, int num_of_labels, 	     
	     __OUT__ SCALAR * means,
	     __OUT__ covs) {
  long i, j;
  int k, m, idx, dim, dim2, *p_str;  
  SCALAR *p_supp, *p_w, *p_mean, *p_cov;
  long *counts;
  
  dim = data->dim;
  dim2 = dim * dim;

  counts = (long*) calloc(num_of_labels, sizeof(int));
  p_supp = data->p_supp; p_w = data->p_w; p_str = data->p_str;

  for (i=0; i<num_of_entries; ++i) {
    _D2_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, 
			 dim, p_str[i], 1., p_supp, dim,
			 p_w, 1,
			 1., means + label[i]*dim, 1);

    p_supp = p_supp + p_str[i] * dim;
    p_w = p_w + p_str[i];
    ++counts[label[i]];
  }

  //  VPRINTF(("Means of different clusters:\n"));
  for (i=0; i<num_of_labels; ++i) {
    _D2_CBLAS_FUNC(scal)(dim, 1./counts[i], means + i*dim, 1);
  }

  p_supp = data->p_supp; p_w = data->p_w; p_str = data->p_str;  
  for (i=0; i<num_of_entries; ++i) {
    p_cov  = covs  + label[i]*dim2;
    p_mean = means + label[i]*dim;
    for (m=0; m<*p_str; ++m) {
      for (j=0,idx=0; j<dim; ++j) for (k=0; k<dim; ++k, ++idx) {
	  p_cov[idx] += *p_w * (p_supp[j]-p_mean[j]) * (p_supp[k]-p_mean[k]);
	}
      p_supp = p_supp + dim;
      p_w = p_w + 1;
    }
    p_str = p_str + 1;    
  }


  //  VPRINTF(("Covariances of different clusters:\n"));
  for (i=0; i<num_of_labels; ++i) {
    _D2_CBLAS_FUNC(scal)(dim2, 1./(counts[i]-dim), covs + i*dim2, 1);
  }
}


  
void d2_mvnrnd(SCALAR * mean, __IN_OUT__ SCALAR * cov, int d, int n, 
	       __OUT__ SCALAR * sample) {
  int info;
  int i,j,k;
  SCALAR *univar_sample;
  
  // Cholesky factorization
  _D2_LAPACKE_FUNC(potrf_)("U", &d, cov, &d, &info); 
  if (info) {VPRINTF(("Cholesky factorization is not successful!")); return;}

  // set lower triangular part to zero
  for (i=0; i<d; ++i) for (j=i+1; j<d; ++j) {
      cov[i*d + j] = 0;
    }

  // Generate std random normal 
  srand (time(NULL));
  univar_sample = _D2_MALLOC_SCALAR(d * n);
  for (i=0; i<d*n; ++i) univar_sample[i] = randn();

  // chol(cov) * univar_sample
  _D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans,  CblasNoTrans,
		       d, n, d, 1., cov, d, univar_sample, d,
		       0.0, sample, d);

  for (i=0,k=0; i<n; ++i) 
    for (j=0; j<d; ++j, ++k) 
      sample[k] += mean[j];
  
}
*/

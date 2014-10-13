//#define _VERBOSE_OUTPUT

#include "d2_math.h"
#include "stdio.h"


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

void d2_mean(sph * data, int * label, int num_of_labels, 	     
	     /** OUT **/ SCALAR * means,
	     /** OUT **/ SCALAR * covs) {
  int i, j, k, m, idx, dim, dim2, *p_str;  
  SCALAR *p_supp, *p_w, *p_mean, *p_cov;
  int *counts;
  
  dim = data->dim;
  dim2 = dim * dim;

  counts = (int*) calloc(num_of_labels, sizeof(int));
  p_supp = data->p_supp; p_w = data->p_w; p_str = data->p_str;
  for (i=0; i<data->size; ++i) {
    
    _D2_CBLAS_FUNC(gemv)(CblasColMajor, CblasNoTrans, 
			 dim, *p_str, 1., p_supp, dim,
			 p_w, 1,
			 1., means + label[i]*dim, 1);

    p_supp = p_supp + *p_str * dim;
    p_w = p_w + *p_str;
    p_str = p_str + 1;
    ++counts[label[i]];
  }

  VPRINTF(("Means of different clusters:\n"));
  for (i=0; i<num_of_labels; ++i) {
    _D2_CBLAS_FUNC(scal)(dim, 1./(counts[i]-dim), means + i*dim, 1);
    VPRINTF(("  mean%d= ", i+1));
    for (j=0; j<dim; ++j) VPRINTF((SCALAR_STDIO_TYPE, means[i*dim+j]));
    VPRINTF(("\n"));
  }

  p_supp = data->p_supp; p_w = data->p_w; p_str = data->p_str;  
  for (i=0; i<data->size; ++i) {
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


  VPRINTF(("Covariances of different clusters:\n"));
  for (i=0; i<num_of_labels; ++i) {
    _D2_CBLAS_FUNC(scal)(dim2, 1./(counts[i]-dim), covs + i*dim2, 1);
    VPRINTF(("  cov%d=", i+1));
    for (j=0; j<dim; ++j) {
      VPRINTF(("\t"));
      for (k=0; k<dim; ++k)
	VPRINTF((SCALAR_STDIO_TYPE, covs[i*dim2+j*dim+k]));
      VPRINTF(("\n"));
    }
  }
}


  
void d2_mvnrnd(SCALAR * mean, /** IN/OUT **/ SCALAR * cov, int d, int n, 
	       /** OUT **/ SCALAR * sample) {
  int info;
  int i,j,k;
  SCALAR *univar_sample;
  
  // Cholesky factorization
  _D2_CLAPACK_FUNC(potrf_)("U", &d, cov, &d, &info); 
  if (info) {VPRINTF(("Cholesky factorization is not successful!")); return;}

  // set lower triangular part to zero
  for (i=0; i<d; ++i) for (j=i+1; j<d; ++j) {
      cov[i*d + j] = 0;
    }

  // Generate std random normal 
  univar_sample = (SCALAR *) malloc(d * n);
  for (i=0; i<d*n; ++i) univar_sample[i] = randn();

  // chol(cov) * univar_sample
  _D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans,  CblasNoTrans,
		       d, n, d, 1., cov, d, univar_sample, d,
		       0.0, sample, d);

  for (i=0,k=0; i<n; ++i) 
    for (j=0; j<d; ++j, ++k) 
      sample[k] += mean[j];
  
}

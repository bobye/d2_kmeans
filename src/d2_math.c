#include "d2_math.h"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif defined __GNUC__
#include <cblas.h>
#include <lapacke.h>
#endif

#ifdef  _D2_DOUBLE
#define _D2_SCALAR          double
#define _D2_FUNC(x)         d ## x
#define _D2_CBLAS_FUNC(x)   cblas_d ## x
#define _D2_LAPACKE_FUNC(x) LAPACKE_d ## x
#elif defined  _D2_SINGLE
#define _D2_SCALAR          float
#define _D2_FUNC(x)         s ## x
#define _D2_CBLAS_FUNC(x)   cblas_s ## x
#define _D2_LAPACKE_FUNC(x) LAPACKE_s ## x
#endif

void d2_mean(sph * data, int * label, int num_of_labels, SCALAR * means) {
  int i, j, dim, *p_str;  
  SCALAR *p_supp, *p_w;
  int *counts;
  
  p_supp = data->p_supp;
  p_w = data->p_w;
  p_str = data->p_str;
  dim = data->dim;

  counts = (int*) calloc(num_of_labels, sizeof(int));

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
    _D2_CBLAS_FUNC(scal)(dim, 1./counts[i], means + i*dim, 1);
    VPRINTF(("\t"));
    for (j=0; j<dim; ++j) printf(SCALAR_STDIO_TYPE, means[i*dim+j]);
    VPRINTF(("\n"));
  }
    

}

void d2_cov(sph * data, int * label, int num_of_labels, SCALAR * means) {
}

void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int n, /** OUT **/ SCALAR * sample) {
}

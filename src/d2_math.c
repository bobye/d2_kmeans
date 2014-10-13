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

void d2_mean(sph * data, int * label, SCALAR * means) {
  int i;

  for (i=0; i<data->size; ++i) {
  }
}

void d2_cov(sph * data, int * label, SCALAR * means) {
}

void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int n, /** OUT **/ SCALAR * sample) {
}

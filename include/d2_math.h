#ifndef _STAT_H_
#define _STAT_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "d2_clustering.h"


#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#elif defined __GNUC__
#include <cblas.h>
#include <clapack.h>
#endif

#ifdef  _D2_DOUBLE
#define _D2_SCALAR          double
#define _D2_FUNC(x)         d ## x
#define _D2_CBLAS_FUNC(x)   cblas_d ## x
#define _D2_CLAPACK_FUNC(x) d ## x
#elif defined  _D2_SINGLE
#define _D2_SCALAR          float
#define _D2_FUNC(x)         s ## x
#define _D2_CBLAS_FUNC(x)   cblas_s ## x
#define _D2_CLAPACK_FUNC(x) s ## x
#endif


  void d2_mean(sph * data, int * label, int num_of_labels, 
	       /** OUT **/ SCALAR * means, /** OUT **/ SCALAR * covs);
  void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int d, int n, /** OUT **/ SCALAR * sample);



#ifdef __cplusplus
}
#endif

#endif /* _STAT_H_ */

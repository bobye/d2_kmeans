#ifndef _STAT_H_
#define _STAT_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "d2_clustering.h"


#include "blas_like.h"



  void d2_mean(sph * data, int * label, int num_of_entries, int num_of_labels, 
	       /** OUT **/ SCALAR * means, /** OUT **/ SCALAR * covs);
  void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int d, int n, /** OUT **/ SCALAR * sample);

  void shuffle(int * array, size_t n);

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) calloc( (x) , sizeof(SCALAR)) 
#define _D2_CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _D2_FREE(x)         free(x)

#elif defined __USE_MKL__
#include <mkl.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) mkl_malloc( (x) *sizeof(SCALAR), 16) 
#define _D2_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) mkl_calloc( (x) , sizeof(SCALAR), 16) 
#define _D2_CALLOC_INT(x)       (int *) mkl_calloc( (x) , sizeof(int), 16)
#define _D2_FREE(x)         mkl_free(x)

#elif defined __GNUC__
#include <cblas.h>
#include <lapacke.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) calloc( (x) , sizeof(SCALAR)) 
#define _D2_CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _D2_FREE(x)         free(x)

#endif


#ifdef __cplusplus
}
#endif


#endif /* _STAT_H_ */

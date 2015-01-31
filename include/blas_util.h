#ifndef _BLAS_UTIL_H_
#define _BLAS_UTIL_H_

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_MALLOC_SIZE_T(x)       (size_t *) malloc( (x) *sizeof(size_t))
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) calloc( (x) , sizeof(SCALAR)) 
#define _D2_CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _D2_CALLOC_SIZE_T(x)       (size_t *) calloc( (x) , sizeof(size_t))
#define _D2_FREE(x)         free(x)

#elif defined __USE_MKL__
#include <mkl.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) mkl_malloc( (x) *sizeof(SCALAR), 16) 
#define _D2_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _D2_MALLOC_SIZE_T(x)       (size_t *) mkl_malloc( (x) *sizeof(size_t), 16)
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) mkl_calloc( (x) , sizeof(SCALAR), 16) 
#define _D2_CALLOC_INT(x)       (int *) mkl_calloc( (x) , sizeof(int), 16)
#define _D2_CALLOC_SIZE_T(x)       (size_t *) mkl_calloc( (x) , sizeof(size_t), 16)
#define _D2_FREE(x)         mkl_free(x)

#elif defined __GNUC__
#include <cblas.h>
#include <lapacke.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_MALLOC_SIZE_T(x)       (size_t *) malloc( (x) *sizeof(size_t))
#define _D2_CALLOC_SCALAR(x)       (SCALAR *) calloc( (x) , sizeof(SCALAR)) 
#define _D2_CALLOC_INT(x)       (int *) calloc( (x) , sizeof(int))
#define _D2_CALLOC_SIZE_T(x)       (size_t *) calloc( (x) , sizeof(size_t))
#define _D2_FREE(x)         free(x)

#endif


#endif /* _BLAS_UTIL_H_ */

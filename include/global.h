#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define _D2_DOUBLE
#define _VERBOSE_OUTPUT

#ifdef _D2_DOUBLE
#define SCALAR double
#define SCALAR_STDIO_TYPE ("%lf ")
#elif defined D2_SINGLE
#define SCALAR float
#define SCALAR_STDIO_TYPE ("%f ")
#endif

#ifdef _VERBOSE_OUTPUT
#define VPRINTF(x) printf x
#else
#define VPRINTF(x) 
#endif

#ifdef  _D2_DOUBLE
#define _D2_SCALAR          double
#define _D2_FUNC(x)         _d ## x
#define _D2_CBLAS_FUNC(x)   cblas_d ## x
#define _D2_CLAPACK_FUNC(x) d ## x
#elif defined  _D2_SINGLE
#define _D2_SCALAR          float
#define _D2_FUNC(x)         _s ## x
#define _D2_CBLAS_FUNC(x)   cblas_s ## x
#define _D2_CLAPACK_FUNC(x) s ## x
#endif

#endif /* _GLOBAL_H_ */

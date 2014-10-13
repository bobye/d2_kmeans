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

#endif /* _GLOBAL_H_ */

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define _D2_DOUBLE

#ifdef _D2_DOUBLE
#define SCALAR double
#define SCALAR_SCANF_TYPE ("%lf")
#elif defined D2_SINGLE
#define SCALAR float
#define SCALAR_SCANF_TYPE ("%f")
#endif

#endif /* _GLOBAL_H_ */

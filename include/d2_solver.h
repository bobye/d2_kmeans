#ifndef _D2_SOLVER_H_
#define _D2_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "blas_like.h"

  void d2_solver_setup();
  void d2_solver_release();
  void d2_solver_debug();

  double d2_match_by_coordinates(int d, int n, double *X, double *wX, int m, double *Y,  double*wY,
				 /** OUT **/ double *x, /** OUT **/ double *lambda);

  double d2_match_by_distmat(int n, int m, double *C, double *wX, double *wY, 
			     /** OUT **/ double *x, /** OUT **/ double *lambda);

#ifdef __cplusplus
}
#endif


#endif /* _D2_SOLVER_H_ */

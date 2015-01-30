#ifndef _D2_SOLVER_H_
#define _D2_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "blas_like.h"

  void d2_solver_setup(size_t num);
  void d2_solver_release();
  void d2_solver_debug();
  /*
  double d2_match_by_coordinates(int d, int n, double *X, double *wX, int m, double *Y,  double*wY,
				 __OUT__ double *x, __OUT__ double *lambda);

  double d2_match_by_symbols(int d, 
			     int n, int *X, double *wX, 
			     int m, int *Y,  double*wY, 
			     int vocab_size, double *dist_mat,
			     __OUT__ double *x, __OUT__ double *lambda);
  */
  double d2_match_by_distmat(int n, int m, double *C, double *wX, double *wY, 
			     /** OUT **/ double *x, /** OUT **/ double *lambda, size_t index);

  double d2_match_by_distmat_qp(int n, int m, double *C, double *L, double rho, double *lw, double *rw, double *x0, /** OUT **/ double *x);
  
  double d2_qpsimple(int str, int count, double *q, /** OUT **/ double *w);
#ifdef __cplusplus
}
#endif


#endif /* _D2_SOLVER_H_ */

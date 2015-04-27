#ifndef _D2_SOLVER_H_
#define _D2_SOLVER_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "common.h"
#include "blas_like.h"

  void d2_solver_setup();
  void d2_solver_release();
  void d2_solver_debug();

  double d2_match_by_distmat(int n, int m, SCALAR *C, SCALAR *wX, SCALAR *wY, 
			     /** OUT **/ SCALAR *x, /** OUT **/ SCALAR *lambda, size_t index);

  double d2_match_by_distmat_qp(int n, int m, SCALAR *C, SCALAR *L, SCALAR rho, SCALAR *lw, SCALAR *rw, SCALAR *x0, /** OUT **/ SCALAR *x);
  
  double d2_qpsimple(int str, int count, SCALAR *q, /** OUT **/ SCALAR *w);
#ifdef __cplusplus
}
#endif


#endif /* _D2_SOLVER_H_ */

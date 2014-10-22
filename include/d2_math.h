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



#ifdef __APPLE__

#include <Accelerate/Accelerate.h>

#elif defined __GNUC__

#include <cblas.h>
#include <clapack.h>

#endif



#ifdef __cplusplus
}
#endif


#endif /* _STAT_H_ */

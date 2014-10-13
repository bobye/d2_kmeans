#ifndef _STAT_H_
#define _STAT_H_


#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "d2_clustering.h"

  void d2_mean(sph * data, int * label, /** OUT **/ SCALAR * means);
  void d2_cov(sph * data, int * label, /** OUT **/ SCALAR * covs);
  void d2_mvnrnd(SCALAR * mean, SCALAR * cov, int n, /** OUT **/ SCALAR * sample);

#ifdef __cplusplus
}
#endif

#endif /* _STAT_H_ */

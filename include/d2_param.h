#ifndef _D2_PARAM_H_
#define _D2_PARAM_H_

typedef struct BADMM_options {
  int maxIters;
  double rhoCoeff;
  int updatePerLoops;
} BADMM_options;

const BADMM_options badmm_clu_options = {500, 2.0, 10};
const BADMM_options badmm_cen_options = {2000, 1.0, 10};

#endif /* _D2_PARAM_H_ */

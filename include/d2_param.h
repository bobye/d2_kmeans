#ifndef _D2_PARAM_H_
#define _D2_PARAM_H_

typedef struct BADMM_options {
  int maxIters;
  double rhoCoeff;
  int updatePerLoops;
} BADMM_options;

#define D2_CENTROID_BADMM    0
#define D2_CENTROID_ADMM     1
#define D2_CENTROID_GRADDEC  2

#endif /* _D2_PARAM_H_ */

#ifndef _D2_PARAM_H_
#define _D2_PARAM_H_

typedef struct {
  int maxIters;
  double rhoCoeff;
  int updatePerLoops;
} BADMM_options;

#define D2_CENTROID_BADMM    (0)
#define D2_CENTROID_ADMM     (1)
#define D2_CENTROID_GRADDEC  (2)


#define D2_EUCLIDEAN_L2      (0)
#define D2_CITYBLOCK_L1      (1)
#define D2_HISTOGRAM         (5)
#define D2_N_GRAM            (6)

#endif /* _D2_PARAM_H_ */

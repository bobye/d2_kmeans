#ifndef _D2_PARAM_H_
#define _D2_PARAM_H_

typedef struct {
  int maxIters;
  double rhoCoeff;
  int updatePerLoops;
} BADMM_options;


typedef struct {
  int maxIters;
  double stepSize;
} GRADDEC_options;


#define D2_CENTROID_BADMM    (0)
#define D2_CENTROID_ADMM     (1)
#define D2_CENTROID_GRADDEC  (2)


/* types of D2 data */
#define D2_EUCLIDEAN_L2      (0)
#define D2_CITYBLOCK_L1      (1)
#define D2_HISTOGRAM         (5)
#define D2_N_GRAM            (6)
#define D2_WORD_EMBED        (7)
#define D2_SPARSE_HISTOGRAM  (12)

#endif /* _D2_PARAM_H_ */

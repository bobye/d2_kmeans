#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"

  typedef struct {
    int dim, str, size;
    int *p_str;
    SCALAR *p_supp;
    SCALAR *p_w;  
  } sph; 

  typedef struct {
    int s_ph, size; // size of phases
    int *label;
    int num_of_labels;
    sph *ph;
  } mph;

  int d2_allocate(mph *p_data,
		  const int size_of_phases,
		  const int size_of_samples,
		  const int *avg_strides,
		  const int *dimension_of_phases);

  int d2_load(void *fp, mph *p_data);

  int d2_free(mph *p_data);


  int d2_centroid_sphBregman(mph *p_data, // data
			     int idx_ph, // index of phases
			     sph *c0,
			     sph *c);

#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

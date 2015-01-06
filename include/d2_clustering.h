#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"

#define __IN__
#define __OUT__ 
#define __IN_OUT__

  // data structure to store d2 of one phase
  typedef struct {
    int dim, str;
    long col, max_col;
    int *p_str;
    long *p_str_cum;
    SCALAR *p_supp;    
    SCALAR *p_w;  
  } sph; 

  // data structure to store d2 of multiple phases
  typedef struct {
    int s_ph /* size of phases */;
    long size /* size of entries */;
    int *label;
    int num_of_labels;
    sph *ph;
  } mph;


  /* basic utilities */
  int d2_allocate_sph(__OUT__ sph *p_data_sph,
		      const int d,
		      const int stride,
		      const long num,
		      const double semicol);

  int d2_allocate(__OUT__ mph *p_data,
		  const int size_of_phases,
		  const long size_of_samples,
		  const int *avg_strides,
		  const int *dimension_of_phases);

  int d2_read(void *fp, __OUT__ mph *p_data);
  int d2_write(void *fp, mph *p_data);
  int d2_free(mph *p_data);

  // working variables that are visible in all algorithms
  typedef struct {
    SCALAR *C;
    SCALAR *X;
    SCALAR *L;
  } var_sph;

  // working variables specific to Bregman ADMM
  typedef struct {
    SCALAR *X, *Z;
    SCALAR *Y;
    SCALAR *Xc, *Zr;
  } var_sphBregman;

  // union of working variables across multiple phases
  typedef struct {
    int s_ph;
    var_sph *g_var;
    var_sphBregman *l_var_sphBregman; // may not initialized, which depends on the actual centroid algorithm used.    
    char *label_switch;
  } var_mph; 


  // interface of random centroids from multivariate normal samples
  int d2_centroid_randn(mph *p_data, 
			int idx_ph, 
			__OUT__ sph *c);

  // interface of random centroids from observations
  int d2_centroid_rands(mph *p_data, 
			int idx_ph, 
			__OUT__ sph *c);
  
  // interface of Bregman ADMM
  int d2_allocate_work_sphBregman(sph *ph, long size,
				  __OUT__ var_sphBregman * var_phwork);
  int d2_free_work_sphBregman(var_sphBregman * var_phwork);
  int d2_centroid_sphBregman(mph *p_data, // data
			     var_mph * var_work, // working data
			     int idx_ph, // index of phases
			     sph *c0,
			     __OUT__ sph *c);


  // interface of Gradient Decent
  int d2_centroid_sphGradDecent(mph *p_data,
				var_mph * var_work,
				int idx_ph,
				sph *c0,
				__OUT__ sph *c);
  
  int d2_centroid_sphADMM(mph *p_data,
			  var_mph *var_work,
			  int idx_ph,
			  sph *c0,
			  __OUT__ sph *c);


  // interface to users
  int d2_clustering(int k, 
		    int max_iter, 
		    mph *p_data, 
		    __OUT__ mph *centroids,
		    int selected_phase);  

#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"

#define __IN__
#define __OUT__ 
#define __IN_OUT__

  /**
   * data structure to store d2 of one phase
   */
  typedef struct {
    int dim, str;
    long col, max_col;
    int *p_str;
    long *p_str_cum;
    SCALAR *p_w;

    /* For data of D2 with Euclidean supports */
    SCALAR *p_supp; 
    
    /**
       For data of D2 with symbolic supports, 
       there is no universal data format specifications.
       One has to write their own IO code for reading/writing such 
       type of data. Here a support is a sequence of fixed length
       @param(dim) and the distance between sequences of symbols
       is defined as the average of their distances along all dimensions,
       where the distance is specified in 2d array @param(dist_mat).
       @param(p_supp_sym) gives the symbolic indices for each dimension. 

       See folder data/dna_seq/ for an example. */
    int *p_supp_sym; 

    /* For data of histograms or D2 with symbolic supports */
    SCALAR *dist_mat; 
  } sph; 


  /**
    data structure for relabeling using triangle inequality:    
    Using the Triangle Inequality to Accelerate k-means, Charles Elkan, ICML03
   */
  typedef struct {
    SCALAR *l; /* lower bound of distance pair */
    SCALAR *u; /* upper bound */
    SCALAR *s; 
    SCALAR *c; /* distance between centroids */
    char *r;    
  } trieq;

  /**
   * data structure to store d2 of multiple phases
   */
  typedef struct {
    int s_ph; /* size of phases */
    long size; /* size of entries */
    int *label;
    int num_of_labels;
    sph *ph;
  } mph;


  /**
   * basic utilities 
   */
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

  int d2_read(const char* filename, __OUT__ mph *p_data);
  int d2_write(const char* filename, mph *p_data);
  int d2_free(mph *p_data);

  /* main algorithm */
  int d2_clustering(int k, 
		    int max_iter, 
		    mph *p_data, 
		    __OUT__ mph *centroids,
		    int selected_phase);  

  /**
   * working variables that are visible in all algorithms
   */
  typedef struct {
    SCALAR *C;
    SCALAR *X;
    SCALAR *L;
  } var_sph;

  /**
   * working variables specific to Bregman ADMM
   */
  typedef struct {
    SCALAR *X, *Z;
    SCALAR *Y;
    SCALAR *Xc, *Zr;
  } var_sphBregman;

  /**
   * union of working variables across multiple phases
   */
  typedef struct {
    int s_ph;
    var_sph *g_var;
    var_sphBregman *l_var_sphBregman; // may not initialized, which depends on the actual centroid algorithm used.    
    char *label_switch;
    trieq tr; /* data structure for relabeling */
  } var_mph; 


  /**
   * interface of random centroids from multivariate normal samples
   */
  int d2_centroid_randn(mph *p_data, 
			int idx_ph, 
			__OUT__ sph *c);

  /**
   * interface of random centroids from observations
   */
  int d2_centroid_rands(mph *p_data, 
			int idx_ph, 
			__OUT__ sph *c);
  
  /**
   * interface of Bregman ADMM
   */
  int d2_allocate_work_sphBregman(sph *ph, long size,
				  __OUT__ var_sphBregman * var_phwork);
  int d2_free_work_sphBregman(var_sphBregman * var_phwork);
  int d2_centroid_sphBregman(mph *p_data, // data
			     var_mph * var_work, // working data
			     int idx_ph, // index of phases
			     sph *c0,
			     __OUT__ sph *c);


  /**
   * interfaces of Gradient Decent and ADMM (for experimental purpose only, deprecated)
   */
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



#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

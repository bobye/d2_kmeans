#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "utils/common.h"
#include "d2/param.h"

  /**
   * data structure to store d2 dataset of one phase
   */
  typedef struct {
    /**
     * This part of data blocks are universal to any data types. 
     * In particular, array @param(p_str, p_str_cum) don't changes after read
     * and @param(p_w) should be kept synchronized.
     *
     * @param(metric_type)
     * D2_EUCLIDEAN_L2 : @param(dim) is the dimension of vector space
     * D2_WORD_EMBED   : @param(dim) is the dimension of embedding space
     * D2_CITYBLOCK_L1 : same above
     * D2_HISTOGRAM    : @param(dim=0)
     * D2_N_GRAM       : @param(dim) is the size of grams.
     * D2_SPARSE_HISTOGRAM : @param(dim=0)
     */
    int dim, str, max_str;
    size_t col, max_col;
    int *p_str;
    size_t *p_str_cum;
    SCALAR *p_w;

    /**
     * Normally, only parts of data blocks are allocated, which is marked as follows. 
     * In MPI implementation, when broadcast data to other nodes, it should firstly
     * check when specific data blocks has been allocated before copying data. 
     * @param(p_supp,p_supp_sym) should be kept synchronized if they are allocated.
     * @param(dist_mat, vocab_size) don't changes after read.
     *
     * what: metric space of supports 
     * @param(metric_type)
     * D2_EUCLIDEAN_L2 : Euclidean space at @param(p_supp,p_w)
     * D2_WORD_EMBED   : Word embedding space (Euclidean)
     * D2_CITYBLOCK_L1 : cityblock metric space @param(p_supp,p_w) => to be implemented
     * D2_HISTOGRAM    : Histogram space at @param(p_w,dist_mat)
     * D2_N_GRAM       : N-gram space at @param(p_supp_sym,dist_mat,vocab_size) 
     *                   and dist_mat is the pairwise distance matrix of vocab
     * D2_SPARSE_HISTOGRAM: Histogram data with sparse non-zero bins 
     *                      @param(p_w,p_sUpp_sym,dist_mat)
     */
    int metric_type;     
    

    /**
     * For data of D2 with Euclidean supports */
    SCALAR *p_supp; 
    
    /**
       Older comments:
       For data of D2 with symbolic supports, 
       there is no universal data format specifications.
       One has to write their own IO code for reading/writing such 
       type of data. Here a support is a sequence of fixed length
       @param(dim) and the distance between sequences of symbols
       is defined as the average of their distances along all dimensions,
       where the distance is specified in 2d array @param(dist_mat).
       @param(p_supp_sym) gives the symbolic indices for each dimension. 

       See folder data/dna_seq/ for an example. 
       
       Comments 2015-02-17
       @param(p_supp_sym) has been used for D2_WORD_EMBED
       
       Comments 2015-07-13
       @param(p_supp_sym) has been used for D2_SPARSE_HISTOGRAM
    */
    int *p_supp_sym; 

    /**
     * Optional @param(dist_mat): 
     * for data in D2_HISTOGRAM, D2_SPARSE_HISTOGRAM and D2_N_GRAM */
    SCALAR *dist_mat; 

    /**
     * Optional @param(vocab_vec): 
     * for data in D2_WORD_EMBED */
    int vocab_size;
    SCALAR *vocab_vec;

    /**
     * Optional @param(is_meta_allocated): 
     * tag to indicate whether vocab_vec or dist_mat is newly allocated */
    char is_meta_allocated;
  } sph; 


  /**
    data structure for relabeling using triangle inequality:    
    Using the Triangle Inequality to Accelerate k-means, Charles Elkan, ICML03
    
    Comments 2015-07-13
    @param(l,u,s) the meaning for each symbol has been consistent with what has
    been described in Elkan's original paper. 
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
    size_t size; /* size of entries */
    size_t global_size;
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
		      const size_t num,
		      const double semicol,
		      const int type);

  int d2_allocate(__OUT__ mph *p_data,
		  const int size_of_phases,
		  const size_t size_of_samples,
		  const int *avg_strides,
		  const int *dimension_of_phases,
		  const int *type_of_phases);

  int d2_read(const char* filename, const char* meta_filename, __OUT__ mph *p_data);
  int d2_write(const char* filename, mph *p_data);
  int d2_write_labels(const char* filename, mph *p_data);
  int d2_write_labels_serial(const char* filename_ind, const char* filename, mph *p_data);

  int d2_write_split(const char* filename, mph *p_data, int splits, char is_pre_processed);
  int d2_free(mph *p_data);

  int d2_init_centroid(mph *p_data, 
		       __OUT__ mph *centroids, 
		       int selected_phase, 
		       int allocate_only);

  /* main algorithm */
  int d2_clustering(int k, 
		    int max_iter, 
		    mph *p_data, 
		    __OUT__ mph *centroids,
		    int selected_phase,
		    char use_triangle,
		    const char* log_file);  

  int d2_assignment(int k,
		    mph *p_data, 
		    mph *centroids, 
		    int selected_phase,
		    const char* centroid_filename,
		    const char* meta_filename);

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

  int d2_allocate_work(mph *p_data, var_mph *var_work, char use_triangle, int selected_phase);
  int d2_free_work(var_mph *var_work, int selected_phase);
  

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
  int d2_allocate_work_sphBregman(sph *ph, size_t size,
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

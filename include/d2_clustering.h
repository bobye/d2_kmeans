#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#define SCALAR double

#ifdef __cplusplus
extern "C" {
#endif

  /* setup problem type of D2 clustering */
  extern int d2_initialize(const int size_of_modalities, 
			   const int *dimensions_of_modalities 
			   );


  /* assign data blocks, memory are pre-allocated */
  extern int d2_assign_data(const int size_of_samples,
			    /** IN **/ int *size_of_supports,
			    /** IN **/ SCALAR *data_block_supp,
			    /** IN **/ SCALAR *data_block_w
			    );

  /* compute the centroid of labeled sub-array, memory for outputs are pre-allocated */
  extern int d2_compute_centroid(const int label_of_interest,
				 /** IN     **/ int *labels,
				 /** IN/OUT **/ int *size_of_supp, 
				 /** OUT    **/ SCALAR *centroid_supp, 
				 /** OUT    **/ SCALAR *centroid_w
				 ); // centroid_supp and centroid_w should allocate enough space subject to size_of_supp

  extern int d2_clustering(const int size_of_clusters,
			   /** OUT **/ int *labels,
			   /** OUT **/ int *size_of_supp,
			   /** OUT **/ SCALAR *cluster_supp,
			   /** OUT **/ SCALAR *cluster_w,
			   /** IN  **/ int reset
			   );

  extern int d2_free();

#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

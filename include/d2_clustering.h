#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#define SCALAR double

#ifdef __cplusplus
extern "C" {
#endif

extern int d2_initialize(const int size_of_modalities,
			 const int *dimensions_of_modalities
			 );



extern int d2_assign_data(const int size_of_samples,
			  /** IN **/ int *size_of_supports,
			  /** IN **/ SCALAR *data_block_supp,
			  /** IN **/ SCALAR *data_block_w
			  );


extern int d2_compute_centroid(const int label_of_interest,
			       /** IN     **/ int *labels,
			       /** IN/OUT **/ int *size_of_supp, 
			       /** OUT    **/ SCALAR *centroid_supp,
			       /** OUT    **/ SCALAR *centroid_w
			       );

extern int d2_clustering(const int size_of_clusters,
			 /** OUT **/ int *labels,
			 /** OUT **/ int *size_of_supp,
			 /** OUT **/ SCALAR *cluster_supp,
			 /** OUT **/ SCALAR *cluster_w
			 );

extern int d2_free();

#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

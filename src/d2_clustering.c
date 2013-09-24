#include <stdlib.h>
#include "d2_clustering.h"


int s_modalities, // size of modalities
  *d_modalities; // dimension of each modality
int s_samples; // size of total samples
int *s_supp; // size of supports, 2D array of s_modalities * s_samples in column major order
SCALAR **p_supp, **p_w; // array of pointers of support vectors and weights, use as p_supp[s_supp[]]


int d2_initialize(const int size_of_modalities,
		  const int *dimensions_of_modalities) {
  int i;
  
  s_modalities = size_of_modalities;

  if (0 >= s_modalities) return -1;// s_modalities should >= 1

  d_modalities = (int*) malloc(s_modalities * sizeof(int));

  if (dimensions_of_modalities == NULL) {
    for (i=0; i<=s_modalities; ++i) d_modalities[i] = 1;
  }
  else {
    for (i=0; i<=s_modalities; ++i) {      
      d_modalities[i] = dimensions_of_modalities[i];
    }
  }

  return 0;
}		  

int d2_free() {
  s_modalities = -1;
  free(d_modalities);
  free(p_supp); free(p_w); 
  return 0;
}

int d2_assign_data(const int size_of_samples,
		   /** IN **/ int *size_of_supports,
		   /** IN **/ SCALAR *data_block_supp,
		   /** IN **/ SCALAR *data_block_w
		   ){
  int i,j,m;
  SCALAR *count_supp = data_block_supp;
  SCALAR *count_w = data_block_w;
    
  s_samples = size_of_samples;
  s_supp = size_of_supports;

  p_supp = ( SCALAR **) malloc(s_samples * s_modalities * sizeof( SCALAR *));
  p_w    = ( SCALAR **) malloc(s_samples * s_modalities * sizeof( SCALAR *));


  j=0; 
  for (i=0; i<s_samples; ++i) 
    for (m=0; m<s_modalities; ++m) {
      p_supp[j] = count_supp;
      p_w[j]    = count_w;
    
      count_supp += s_supp[j]* d_modalities[m];
      count_w += s_supp[j];
      ++j;
    }

  return 0;

}

int d2_compute_centroid(const int label_of_interest,
			/** IN     **/ int *labels,
			/** IN/OUT **/ int *size_of_supp, 
			/** OUT    **/ SCALAR *centroid_supp, 
			/** OUT    **/ SCALAR *centroid_w,
			/** IN     **/ int reset
			){
  // step 1: setup initial guess
  if (reset == 1) {
    int i, nzm =0, nzc =0;
    SCALAR *pca_mat, *centers;
    for (int i=0; i<s_modalities; ++i) {
      nzc += d_modalities[i];
      nzm += d_modalities[i] * d_modalities[i];
    }
    pca_mat = (SCALAR *) malloc( nzm * sizeof(SCALAR) );
    // PCA setup
    for (i=0; i<size_of_samples; ++i)
      if (-1 == label_of_interest || labels[i] == label_of_interest) {
	
      } 
  }
  // step 2a: fix centroid_supp, update centroid_w
  // step 2b: fix centroid_w, update centroid_supp
  // step 3: stop criterion, output statistics
}


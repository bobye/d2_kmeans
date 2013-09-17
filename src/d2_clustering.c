#include <stdlib.h>
#include "d2_clustering.h"


int s_modalities, *d_modalities;
int s_samples;
int *s_supp;
SCALAR **p_supp, **p_w;

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
		   int *size_of_supports,
		   SCALAR *data_block_supp,
		   SCALAR *data_block_w
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
      ++j;
      p_supp[j] = count_supp;
      p_w[j]    = count_w;
    
      count_supp += s_supp[j]* d_modalities[m];
      count_w += s_supp[j];
    }

  return 0;

}


#include "d2_clustering.h"
#include "d2_math.h"

int d2_centroid_randn(mph *p_data, int idx_ph, sph *c) {
  int i;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  
  // initialization 
  d2_allocate_sph(c, dim, str, num_of_labels, 0.);
  
  // set stride
  for (i=0; i<num_of_labels; ++i) {
    c->p_str[i] = str;
  }
  
  // set weight
  for (i=0; i<num_of_labels*str; ++i) {
    c->p_w[i] = 1./str;
  }

  // set column
  c->col = str * num_of_labels;
  
  // compute mean and cov
  SCALAR *means, *covs;
  means = (SCALAR *) calloc(dim * num_of_labels, sizeof(SCALAR));
  covs  = (SCALAR *) calloc(dim * dim * num_of_labels, sizeof(SCALAR));
  d2_mean(data_ph, p_data->label, p_data->size, num_of_labels, means, covs);   
  
  // generate multivariate normal
  for (i=0; i<num_of_labels; ++i) {      
    d2_mvnrnd(means+i*dim, covs+i*dim*dim, dim, str, c->p_supp + i*str*dim);
  }
  free(means);
  free(covs);
  return 0;
}

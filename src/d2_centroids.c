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
  
  // compute mean and cov
  SCALAR *means, *covs;
  means = (SCALAR *) calloc(dim * num_of_labels, sizeof(SCALAR));
  covs  = (SCALAR *) calloc(dim * dim * num_of_labels, sizeof(SCALAR));
  d2_mean(data_ph, p_data->label, num_of_labels, means, covs);   
  
  // generate multivariate normal
  for (i=0; i<num_of_labels; ++i) {      
    d2_mvnrnd(means+i*dim, covs+i*dim*dim, dim, str, c->p_supp + i*str*dim);
  }
}

int d2_centroid_sphBregman(mph *p_data, // data
			   int idx_ph, // index of phases
			   sph *c0,
			   /** OUT **/ sph *c) {
  int i,j,k;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  int *p_str = data_ph->p_str;
  SCALAR *p_supp = data_ph->p_supp;
  SCALAR *p_w = data_ph->p_w;
  
  if (!c0) {
    d2_centroid_randn(p_data, idx_ph, c);
  } else {
    *c = *c0; // warm start
  }

  



  return 0;
}





int d2_centroid_sphGD     (mph *p_data,
			   int idx_ph,
			   sph *c0,
			   /** OUT **/ sph *c) {
  return 0;
}


int d2_centroid_sphADMM    (mph *p_data,
			    int idx_ph,
			    sph *c0,
			    /** OUT **/ sph *c) {
  return 0;
}

int d2_centroid_sphLP      (mph *p_data,
			    int idx_ph,
			    sph *c0,
			    /** OUT **/ sph *c) {
  return 0;
}

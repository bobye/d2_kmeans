#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"
#include "d2_math.h"

int d2_free(mph *p_data) {
  int i;
  for (i=0; i<p_data->s_ph; ++i) {
    free(p_data->ph[i].p_str);
    free(p_data->ph[i].p_supp);
    free(p_data->ph[i].p_w);
  }
  free(p_data->ph);
  return 0;
}

int d2_allocate_sph(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const int num,
		    const double semicol) {

  int n;
  n = num * d * (stride + semicol);

  p_data_sph->dim = d;
  p_data_sph->str = stride;
  p_data_sph->size = num;

  p_data_sph->p_str  = (int *) malloc(num * sizeof(int));
  p_data_sph->p_supp = (SCALAR *) malloc(n * sizeof(SCALAR));
  p_data_sph->p_w    = (SCALAR *) malloc(n * sizeof(SCALAR));  

  if (!(p_data_sph->p_str && p_data_sph->p_supp && p_data_sph->p_w)) {
    return -1;
  }

  return 0;
}

int d2_free_sph(sph *p_data_sph) {
  free(p_data_sph->p_str);
  free(p_data_sph->p_supp);
  free(p_data_sph->p_w);
  return 0;
}

int d2_allocate(mph *p_data,
		const int size_of_phases,
		const int size_of_samples,
		const int *avg_strides,
		const int *dimension_of_phases) {
  int i;
  int success = 0;

  p_data->s_ph = size_of_phases;
  p_data->size = size_of_samples; 
  p_data->ph   = (sph *) malloc(size_of_phases * sizeof(sph));


  for (i=0; i<p_data->s_ph; ++i) {

    p_data->label = (int *) calloc(size_of_samples, sizeof(int)); // initialize to zero
    success = d2_allocate_sph(p_data->ph + i, 
			      dimension_of_phases[i], 
			      avg_strides[i], 
			      size_of_samples, 
			      0.6);
    if (success != 0) break;
  }
  return success;
}

int d2_load(void *fp_void, mph *p_data) {
  FILE *fp = (FILE *) fp_void;

  int i,j,n,dim;
  int **p_str, str, strxdim;
  double **p_supp, **p_w, *p_supp_sph, *p_w_sph;
  int s_ph = p_data->s_ph;
  int size = p_data->size;

  p_str  = (int **) malloc(s_ph * sizeof(int *));
  p_supp = (double **) malloc(s_ph * sizeof(double *));
  p_w    = (double **) malloc(s_ph * sizeof(double *));

  for (n=0; n<s_ph; ++n) {
    p_str[n]  = p_data->ph[n].p_str;
    p_supp[n] = p_data->ph[n].p_supp;
    p_w[n]    = p_data->ph[n].p_w;
  }

  for (i=0; i<size; ++i) {
    for (n=0; n<s_ph; ++n) {      
      // read dimension and stride    
      fscanf(fp, "%d", &dim); assert(dim == p_data->ph[n].dim);
      fscanf(fp, "%d", p_str[n]); 
      str = *(p_str[n]);

      // read weights      
      p_supp_sph = p_supp[n];
      for (j=0; j<str; ++j, ++p_supp_sph)
	fscanf(fp, SCALAR_SCANF_TYPE, p_supp_sph); 
      p_supp[n] = p_supp_sph;

      // read support vec
      p_w_sph = p_w[n];strxdim = str*dim;
      for (j=0; j<strxdim; ++j, ++p_w_sph)
	fscanf(fp, SCALAR_SCANF_TYPE, p_w_sph); 
      p_w[n] = p_w_sph;

      p_str[n] ++;
    }
  }

  return 0;
}



int d2_centroid_sphBregman(mph *p_data, // data
			   int idx_ph, // index of phases
			   sph *c0) {
  int i,j,k;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  int *p_str = data_ph->p_str;
  SCALAR *p_supp = data_ph->p_supp;
  SCALAR *p_w = data_ph->p_w;


  if (!c0) {
    // initialization 
    d2_allocate_sph(c0, dim, str, num_of_labels, 0);
    for (i=0; i<num_of_labels; ++i) {
      c0->p_str[i] = str;
    }
    for (i=0; i<num_of_labels*str; ++i) {
      c0->p_w[i] = 1./str;
    }

    // compute mean and cov
    SCALAR *means, *covs;
    means = (SCALAR *) malloc(dim * num_of_labels);
    covs  = (SCALAR *) malloc(dim * dim * num_of_labels);
    d2_mean(data_ph, p_data->label, means);
    d2_cov(data_ph, p_data->label, covs);
    for (i=0; i<num_of_labels; ++i) {      
      d2_mvnrnd(means+i*dim, covs+i*dim*dim, str, c0->p_supp + i*str*dim);
    }
  }
}


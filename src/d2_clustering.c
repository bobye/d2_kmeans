#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"


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

int d2_allocate(mph *p_data,
		const int size_of_phases,
		const int size_of_samples,
		const int *avg_strides,
		const int *dimension_of_phases) {
  int i,n;
  int success;

  p_data->s_ph = size_of_phases;
  p_data->size = size_of_samples; 
  p_data->ph   = (sph *) malloc(size_of_phases * sizeof(sph));

  success = 0;
  for (i=0; i<p_data->s_ph; ++i) {

    p_data->ph[i].dim     = dimension_of_phases[i];
    p_data->ph[i].avg_str = avg_strides[i];

    n = p_data->size * p_data->ph[i].dim * (p_data->ph[i].avg_str + 0.6);

    p_data->ph[i].p_str  = (int *) malloc(size_of_samples * sizeof(int)); 
    p_data->ph[i].p_supp = (SCALAR *) malloc(n * sizeof(SCALAR));
    p_data->ph[i].p_w    = (SCALAR *) malloc(n * sizeof(SCALAR));

    if (!(p_data->ph[i].p_str && p_data->ph[i].p_supp && p_data->ph[i].p_w)) {
      success=-1;
      break;
    }
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




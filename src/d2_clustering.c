#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"
#include "d2_math.h"


extern int d2_alg_type;

int d2_free(mph *p_data) {
  int i;
  for (i=0; i<p_data->s_ph; ++i) {
    _D2_FREE(p_data->ph[i].p_str);
    _D2_FREE(p_data->ph[i].p_supp);
    _D2_FREE(p_data->ph[i].p_w);
  }
  _D2_FREE(p_data->ph);
  if (!p_data->label) _D2_FREE(p_data->label);
  return 0;
}

int d2_allocate_sph(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const int num,
		    const double semicol) {

  int n, m;
  assert(d>0 && stride >0 && num >0);

  n = num * (stride + semicol) * d; // pre-allocate
  m = num * (stride + semicol); p_data_sph->max_col = m;

  p_data_sph->dim = d;  
  p_data_sph->str = stride;
  //  p_data_sph->size = num;


  p_data_sph->p_str  = _D2_CALLOC_INT(num);
  p_data_sph->p_supp = _D2_MALLOC_SCALAR(n);
  p_data_sph->p_w    = _D2_MALLOC_SCALAR(m);


  if (!(p_data_sph->p_str && p_data_sph->p_supp && p_data_sph->p_w)) {
    return -1;
  }

  return 0;
}

int d2_free_sph(sph *p_data_sph) {
  _D2_FREE(p_data_sph->p_str);
  _D2_FREE(p_data_sph->p_supp);
  _D2_FREE(p_data_sph->p_w);
  return 0;
}

int d2_allocate(mph *p_data,
		const int size_of_phases,
		const int size_of_samples,
		const int *avg_strides, /* it is very important to make sure 
					   that avg_strides are specified correctly,
					   otherwise it will cause memory problems. */
		const int *dimension_of_phases) {
  int i;
  int success = 0;

  p_data->s_ph = size_of_phases;
  p_data->size = size_of_samples; 
  p_data->ph   = (sph *) malloc(size_of_phases * sizeof(sph));


  for (i=0; i<p_data->s_ph; ++i) {

    p_data->label = _D2_CALLOC_INT(size_of_samples); // initialize to zero
    p_data->num_of_labels = 1; // default
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

  int i,j,n,dim,c;
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
      c=fscanf(fp, "%d", &dim); 
      if (c!=1) {
	VPRINTF(("Warning: only read %d d2!\n", i));
	p_data->size = i;
	_D2_FREE(p_w); _D2_FREE(p_supp); _D2_FREE(p_str);
	return 0;
      }
      assert(dim == p_data->ph[n].dim);
      fscanf(fp, "%d", p_str[n]); 
      str = *(p_str[n]); assert(str > 0);
      p_data->ph[n].col += str; assert(p_data->ph[n].col < p_data->ph[n].max_col);

      // read weights      
      p_w_sph = p_w[n];
      for (j=0; j<str; ++j)
	fscanf(fp, SCALAR_STDIO_TYPE, &p_w_sph[j]); 
      p_w[n] = p_w[n] + str;

      // read support vec
      p_supp_sph = p_supp[n];strxdim = str*dim;
      for (j=0; j<strxdim; ++j)
	fscanf(fp, SCALAR_STDIO_TYPE, &p_supp_sph[j]); 
      p_supp[n] = p_supp[n] + strxdim;

      p_str[n] ++;
    }
  }

  _D2_FREE(p_w); _D2_FREE(p_supp); _D2_FREE(p_str);

  return 0;
}

int d2_allocate_work(mph *p_data, var_mph *var_work) {
  int i;
  var_work->s_ph = p_data->s_ph;

  var_work->g_var = (var_sph *) malloc(p_data->s_ph * sizeof(var_sph));
  if (d2_alg_type == 0) 
    var_work->l_var_sphBregman = (var_sphBregman *) malloc(p_data->s_ph * sizeof(var_sphBregman));

  for (i=0; i<p_data->s_ph; ++i) {
    var_work->g_var[i].C = 
      _D2_MALLOC_SCALAR(p_data->ph[i].str * p_data->ph[i].col);

    if (d2_alg_type == 0) {
      d2_allocate_work_sphBregman(p_data->ph +i, p_data->size, 
				  var_work->l_var_sphBregman+i);
    }
  }
  return 0;
}

int d2_free_work(var_mph *var_work) {
  int i;
  for (i=0; i<var_work->s_ph; ++i) {
    _D2_FREE(var_work->g_var[i].C);
    if (d2_alg_type) {
      d2_free_work_sphBregman(var_work->l_var_sphBregman + i);
    }
  }
  free(var_work->g_var);
  free(var_work->l_var_sphBregman);
  return 0;
}

int d2_labeling(mph *p_data,
		mph *centroids,
		int num_clusters) {

  return 0;
}

int d2_clustering(int num_clusters, 
		  int max_iter, 
		  mph *p_data, 
		  /** OUT **/ mph *centroids){
  int i,j,k,iter;
  int s_ph = p_data->s_ph;
  int size = p_data->size;
  int *label = p_data->label;
  assert(num_clusters>0 && max_iter > 0);

  for (i=0; i<size; ++i) 
    label[i] = rand() % num_clusters;

  // initialize centroids from random
  centroids->s_ph = s_ph;
  centroids->size = num_clusters;
  centroids->ph = (sph *) malloc(s_ph * sizeof(sph));
  for (i=0; i<s_ph; ++i)
    d2_centroid_randn(p_data, i, centroids->ph + i);

  // initialize auxiliary variables
  var_mph var_work;
  d2_allocate_work(p_data, &var_work);

  // start centroid-based clustering here
  for (iter=0; iter<max_iter; ++iter) {
    VPRINTF(("Round %d ... \n", iter));
    d2_labeling(p_data, centroids, num_clusters);

    for (i=0; i<s_ph; ++i) {
      VPRINTF(("\t phase %d: \n", i));            

      if (d2_alg_type == 0) 
	d2_centroid_sphBregman(p_data, &var_work, i, centroids->ph + i, centroids->ph + i);

    }
  }

  d2_free_work(&var_work);

  return 0;
}

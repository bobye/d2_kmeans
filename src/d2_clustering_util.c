#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include "d2_clustering.h"
#include "d2_math.h"
#include "d2_param.h"

extern int d2_alg_type;


/**
 * Allocate memory for a single phase: 
 * @param(d): if (d == 0) then only allocate space for weights
 * @param(semicol): auxiliary variable that request extra spaces, 0<=semicol<1
 */
int d2_allocate_sph(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const size_t num,
		    double semicol) {

  size_t n, m;
  assert(stride >0 && num >0 && semicol >= 0);

  n = num * (stride + semicol) * d; // pre-allocate
  m = num * (stride + semicol); p_data_sph->max_col = m;

  p_data_sph->dim = d;  
  p_data_sph->str = stride;
  //  p_data_sph->size = num;


  p_data_sph->p_str  = _D2_CALLOC_INT(num);
  p_data_sph->p_str_cum  = _D2_CALLOC_SIZE_T(num);
  p_data_sph->p_w    = _D2_MALLOC_SCALAR(m);

  // consider different data format
  if (d>0) {
    p_data_sph->p_supp = _D2_MALLOC_SCALAR(n);
    p_data_sph->metric_type = D2_EUCLIDEAN_L2;
  }
  else if (d==0) { 
    p_data_sph->dist_mat = _D2_MALLOC_SCALAR( stride * stride );
    p_data_sph->metric_type = D2_HISTOGRAM;
  }

  p_data_sph->p_supp_sym = NULL;

  return 0;
}

/**
 * Free the space of a single phase
 */
int d2_free_sph(sph *p_data_sph) {
  _D2_FREE(p_data_sph->p_w);
  _D2_FREE(p_data_sph->p_str);
  _D2_FREE(p_data_sph->p_str_cum);

  if (p_data_sph->metric_type == D2_EUCLIDEAN_L2 ||
      p_data_sph->metric_type == D2_CITYBLOCK_L1) {
    _D2_FREE(p_data_sph->p_supp);
  }
  else if (p_data_sph->metric_type == D2_HISTOGRAM) {
    _D2_FREE(p_data_sph->dist_mat);
  }
  else if (p_data_sph->metric_type == D2_N_GRAM) {
    _D2_FREE(p_data_sph->p_supp_sym);
    _D2_FREE(p_data_sph->dist_mat);
  }
  return 0;
}

/**
 * Allocate memory for data, it is possible that the pre-allocated memory is
 * insufficient when loading data. In that case, memory will be reallocated. 
 */
int d2_allocate(mph *p_data,
		const int size_of_phases,
		const size_t size_of_samples,
		const int *avg_strides, /**
					   It is very important to make sure 
					   that avg_strides are specified correctly.
					   It articulates how sparse the centroid could be. 
					   By default, it should be the average number
					   of bins of data objects.
					*/
		const int *dimension_of_phases) {
  size_t i;
  int success = 0;

  p_data->s_ph = size_of_phases;
  p_data->size = size_of_samples; 
  p_data->ph   = (sph *) malloc(size_of_phases * sizeof(sph));
  p_data->num_of_labels = 0; // default

  // initialize to all labels to invalid -1
  p_data->label = _D2_MALLOC_INT(size_of_samples); 
  for (i=0; i<p_data->size; ++i)  p_data->label[i] = -1;

  for (i=0; i<p_data->s_ph; ++i) {
    success = d2_allocate_sph(p_data->ph + i, 
			      dimension_of_phases[i], 
			      avg_strides[i], 
			      size_of_samples, 
			      0.6);
    if (success != 0) break;
  }

  return success;
}

/**
 * Free the space of entire data
 */
int d2_free(mph *p_data) {
  int i;
  for (i=0; i<p_data->s_ph; ++i) {
    if (p_data->ph[i].col > 0) d2_free_sph(p_data->ph + i);
  }
  free(p_data->ph);
  if (!p_data->label) _D2_FREE(p_data->label);
  return 0;
}

/**
 * Allocate memory for working data
 */
int d2_allocate_work(mph *p_data, var_mph *var_work, char use_triangle) {
  int i;
  size_t size = p_data->size;
  int num_of_labels = p_data->num_of_labels;
  trieq *p_tr = &var_work->tr;
  var_work->s_ph = p_data->s_ph;

  var_work->g_var = (var_sph *) malloc(p_data->s_ph * sizeof(var_sph));
  if (d2_alg_type == D2_CENTROID_BADMM) {
      var_work->l_var_sphBregman = (var_sphBregman *) malloc(p_data->s_ph * sizeof(var_sphBregman));
      assert(var_work->l_var_sphBregman);
  }
  for (i=0; i<p_data->s_ph; ++i) {
    int str = p_data->ph[i].str;
    int col = p_data->ph[i].col;

    var_work->g_var[i].C = NULL;
    var_work->g_var[i].X = NULL;
    var_work->g_var[i].L = NULL;

    var_work->g_var[i].C = _D2_MALLOC_SCALAR(str * (col + num_of_labels*str)); // space for transportation cost
    assert(var_work->g_var[i].C);

    if (d2_alg_type == D2_CENTROID_BADMM) {
      d2_allocate_work_sphBregman(p_data->ph +i, p_data->size, 
				  var_work->l_var_sphBregman+i);
    }
    if (d2_alg_type == D2_CENTROID_ADMM) {
      var_work->g_var[i].X = _D2_MALLOC_SCALAR(str * col);
      assert(var_work->g_var[i].X);
    }
    if (d2_alg_type == D2_CENTROID_GRADDEC) {
      var_work->g_var[i].X = _D2_MALLOC_SCALAR(str * (col + num_of_labels * str));
      var_work->g_var[i].L = _D2_MALLOC_SCALAR(str * (size + num_of_labels));
      assert(var_work->g_var[i].X);
      assert(var_work->g_var[i].L);
    }
  }
  var_work->label_switch = (char *) malloc(size * sizeof(char)); 

  if (use_triangle) {
    size_t j;
    p_tr->l = _D2_MALLOC_SCALAR(size * num_of_labels);
    p_tr->u = _D2_MALLOC_SCALAR(size);
    p_tr->s = _D2_MALLOC_SCALAR(num_of_labels);
    p_tr->c = _D2_MALLOC_SCALAR(num_of_labels * num_of_labels);
    p_tr->r = (char *) calloc(size, sizeof(char));

    for (j=0; j<size * num_of_labels; ++j) p_tr->l[j] = 0;
    for (j=0; j<size; ++j) {p_tr->u[j] = DBL_MAX; p_tr->r[j] = 1; }
  }
  return 0;
}

/**
 * Free space for working data
 */
int d2_free_work(var_mph *var_work) {
  int i;
  trieq *p_tr = &var_work->tr;

  for (i=0; i<var_work->s_ph; ++i) {
    if (var_work->g_var[i].C) _D2_FREE(var_work->g_var[i].C);
    if (var_work->g_var[i].X) _D2_FREE(var_work->g_var[i].X);
    if (var_work->g_var[i].L) _D2_FREE(var_work->g_var[i].L);

    if (d2_alg_type == D2_CENTROID_BADMM) {
      d2_free_work_sphBregman(var_work->l_var_sphBregman + i);
    }
  }

  if (var_work->g_var) free(var_work->g_var);
  if (d2_alg_type == D2_CENTROID_BADMM) free(var_work->l_var_sphBregman);
  if (var_work->label_switch) free(var_work->label_switch);
  if (p_tr->l) _D2_FREE(p_tr->l);
  if (p_tr->u) _D2_FREE(p_tr->u);
  if (p_tr->s) _D2_FREE(p_tr->s);
  if (p_tr->c) _D2_FREE(p_tr->c);
  if (p_tr->r) free(p_tr->r);
  return 0;
}


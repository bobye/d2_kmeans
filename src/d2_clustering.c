#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"
#include "d2_math.h"
#include "d2_solver.h"
#include "d2_param.h"

#include <omp.h>

extern int d2_alg_type;

int d2_free(mph *p_data) {
  int i;
  for (i=0; i<p_data->s_ph; ++i) {
    _D2_FREE(p_data->ph[i].p_str);
    _D2_FREE(p_data->ph[i].p_supp);
    _D2_FREE(p_data->ph[i].p_w);
  }
  free(p_data->ph);
  if (!p_data->label) _D2_FREE(p_data->label);
  return 0;
}

int d2_allocate_sph(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const long num,
		    const double semicol) {

  int n, m;
  assert(d>0 && stride >0 && num >0);

  n = num * (stride + semicol) * d; // pre-allocate
  m = num * (stride + semicol); p_data_sph->max_col = m;

  p_data_sph->dim = d;  
  p_data_sph->str = stride;
  //  p_data_sph->size = num;


  p_data_sph->p_str  = _D2_CALLOC_INT(num);
  p_data_sph->p_str_cum  = _D2_CALLOC_LONG(num);
  p_data_sph->p_supp = _D2_MALLOC_SCALAR(n);
  p_data_sph->p_w    = _D2_MALLOC_SCALAR(m);


  if (!(p_data_sph->p_str && p_data_sph->p_supp && p_data_sph->p_w)) {
    return -1;
  }

  return 0;
}

int d2_free_sph(sph *p_data_sph) {
  _D2_FREE(p_data_sph->p_str);
  _D2_FREE(p_data_sph->p_str_cum);
  _D2_FREE(p_data_sph->p_supp);
  _D2_FREE(p_data_sph->p_w);
  return 0;
}

/** Allocate memory for data, it is possible that the pre-allocated memory is
    insufficient when loading data. In that case, memory will be reallocated. 
 */
int d2_allocate(mph *p_data,
		const int size_of_phases,
		const long size_of_samples,
		const int *avg_strides, /* It is very important to make sure 
					   that avg_strides are specified correctly.
					   It articulates how sparse the centroid could be. 
					   By default, it should be the average number
					   of bins of data objects.
					*/
		const int *dimension_of_phases) {
  int i;
  int success = 0;

  p_data->s_ph = size_of_phases;
  p_data->size = size_of_samples; 
  p_data->ph   = (sph *) malloc(size_of_phases * sizeof(sph));
  p_data->num_of_labels = 0; // default

  // initialize to all labels to invalid -1
  p_data->label = _D2_MALLOC_INT(size_of_samples); 
  //  for (i=0; i<p_data->size; ++i)  p_data->label[i] = -1;

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


int d2_allocate_work(mph *p_data, var_mph *var_work) {
  int i;
  var_work->s_ph = p_data->s_ph;

  var_work->g_var = (var_sph *) malloc(p_data->s_ph * sizeof(var_sph));
  if (d2_alg_type == D2_CENTROID_BADMM) 
    var_work->l_var_sphBregman = (var_sphBregman *) malloc(p_data->s_ph * sizeof(var_sphBregman));

  for (i=0; i<p_data->s_ph; ++i) {
    var_work->g_var[i].C = 
      _D2_MALLOC_SCALAR(p_data->ph[i].str * p_data->ph[i].col);

    var_work->g_var[i].X = NULL;
    var_work->g_var[i].L = NULL;

    if (d2_alg_type == D2_CENTROID_BADMM) {
      d2_allocate_work_sphBregman(p_data->ph +i, p_data->size, 
				  var_work->l_var_sphBregman+i);
    }
    if (d2_alg_type == D2_CENTROID_ADMM) {
      var_work->g_var[i].X = 
	_D2_MALLOC_SCALAR(p_data->ph[i].str * p_data->ph[i].col);
    }
    if (d2_alg_type == D2_CENTROID_GRADDEC) {
      var_work->g_var[i].X = 
	_D2_MALLOC_SCALAR(p_data->ph[i].str * p_data->ph[i].col);
      var_work->g_var[i].L = 
	_D2_MALLOC_SCALAR(p_data->ph[i].str * p_data->size);
    }
  }

  var_work->label_switch = (char *) malloc(p_data->size * sizeof(char)); 

  return 0;
}

int d2_free_work(var_mph *var_work) {
  int i;
  for (i=0; i<var_work->s_ph; ++i) {
    _D2_FREE(var_work->g_var[i].C);
    if (var_work->g_var[i].X) _D2_FREE(var_work->g_var[i].X);
    if (var_work->g_var[i].L) _D2_FREE(var_work->g_var[i].L);

    if (d2_alg_type == D2_CENTROID_BADMM) {
      d2_free_work_sphBregman(var_work->l_var_sphBregman + i);
    }
  }
  free(var_work->g_var);
  free(var_work->l_var_sphBregman);
  free(var_work->label_switch);
  return 0;
}

/** Compute the distance from each point to the all centroids.
    This part can be parallelized.
 */
long d2_labeling(__IN_OUT__ mph *p_data,
		mph *centroids,
		var_mph * var_work,
		bool isFirstTime,
		int selected_phase) {
  long i, count = 0;
  int **p_str;
  int n;
  double **p_supp, **p_w;
  long **p_str_cum;
  int s_ph = p_data->s_ph;
  long size = p_data->size;

  p_str  = (int **) malloc(s_ph * sizeof(int *));
  p_supp = (double **) malloc(s_ph * sizeof(double *));
  p_w    = (double **) malloc(s_ph * sizeof(double *));
  p_str_cum = (long **) malloc(s_ph * sizeof(long *));

  for (n=0; n<s_ph; ++n) {
    p_str[n]  = p_data->ph[n].p_str;
    p_supp[n] = p_data->ph[n].p_supp;
    p_w[n]    = p_data->ph[n].p_w;
    p_str_cum[n] = p_data->ph[n].p_str_cum;
  }

#pragma omp parallel for  
  for (i=0; i<size; ++i) {
    double min_distance = -1;	
    int min_label = -1;
    int n;
    long j;

    for (j=0; j<centroids->size; ++j) {
      double d = 0.0, val;
      for (n=0; n<p_data->s_ph; ++n) 
	if (selected_phase < 0 || n == selected_phase) {
	  int str = centroids->ph[n].str;
	  int dim = p_data->ph[n].dim;
	  assert(dim == centroids->ph[n].dim);
	  if (d2_alg_type == D2_CENTROID_BADMM) { 
	    val = d2_match_by_coordinates(dim, 
					  p_str[n][i], p_supp[n] + p_str_cum[n][i]*dim, p_w[n] + p_str_cum[n][i],
					  str, centroids->ph[n].p_supp + j*str*dim, centroids->ph[n].p_w + j*str, 
					  NULL, // x and lambda are implemented later
					  NULL);
	    d += val;
	  }
	  if (d2_alg_type == D2_CENTROID_GRADDEC) {
	    val = d2_match_by_coordinates(dim, 					  
					  str, centroids->ph[n].p_supp + j*str*dim, centroids->ph[n].p_w + j*str, 
					  p_str[n][i], p_supp[n] + p_str_cum[n][i]*dim, p_w[n] + p_str_cum[n][i],
					  var_work->g_var[n].X + p_str_cum[n][i]*str, 
					  var_work->g_var[n].L + i*str);
	  }
	}

      if (min_distance < 0 || d < min_distance) {
	min_distance = d; min_label = j;
      }
    }

    if (p_data->label[i] == min_label && !isFirstTime) {
      if (d2_alg_type == 0) {
	var_work->label_switch[i] = 0;
      }
    } else {
      p_data->label[i] = min_label;
      if (d2_alg_type == 0) {
	var_work->label_switch[i] = 1;
      }
      count ++;
    }
  }

  free(p_str); free(p_supp); free(p_w); free(p_str_cum);
  VPRINTF(("\t %ld objects change their labels [done]\n", count));
  
  return count;
}


/** the main algorithm for d2 clustering */
int d2_clustering(int num_of_clusters, 
		  int max_iter, 
		  mph *p_data, 
		  __OUT__ mph *centroids,
		  int selected_phase){
  long i;
  int iter;
  int s_ph = p_data->s_ph;
  long size = p_data->size;
  int *label = p_data->label;
  long label_change_count;
  var_mph var_work;

  assert(num_of_clusters>0 && max_iter > 0);

  // label all objects as invalid numbers
  p_data->num_of_labels = num_of_clusters;
  for (i=0; i<size; ++i) label[i] = rand() % num_of_clusters;

  // initialize centroids from random
  centroids->s_ph = s_ph;
  centroids->size = num_of_clusters;
  centroids->ph = (sph *) malloc(s_ph * sizeof(sph));
  for (i=0; i<s_ph; ++i) 
    if (selected_phase < 0 || i == selected_phase) 
      d2_centroid_rands(p_data, i, centroids->ph + i);
  // d2_write(stdout, centroids); getchar();

  // initialize auxiliary variables
  d2_allocate_work(p_data, &var_work);

  // start centroid-based clustering here
  for (iter=0; iter<max_iter; ++iter) {
    VPRINTF(("Round %d ... \n", iter));
    VPRINTF(("\tRe-labeling all instances ... ")); VFLUSH;
    d2_solver_setup();
    if (iter == 0)
      label_change_count = d2_labeling(p_data, centroids, &var_work, true, selected_phase);
    else 
      label_change_count = d2_labeling(p_data, centroids, &var_work, false, selected_phase);
    d2_solver_release();

    /* termination criterion */
    if (label_change_count < 0.005 * size) break;

    VPRINTF(("\tUpdate centroids ... \n"));

    for (i=0; i<s_ph; ++i) 
      if (selected_phase < 0 || i == selected_phase) {
	VPRINTF(("\t phase %ld: \n", i));            
      
	if (d2_alg_type == 0) 
	  d2_centroid_sphBregman(p_data, &var_work, i, centroids->ph + i, centroids->ph + i);
      }
  }

  d2_free_work(&var_work);

  return 0;
}


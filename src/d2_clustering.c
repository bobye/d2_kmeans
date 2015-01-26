#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include "d2_clustering.h"
#include "d2_math.h"
#include "d2_solver.h"
#include "d2_param.h"

#ifndef __APPLE__
#include <omp.h>
#endif 

extern int d2_alg_type;

int d2_allocate_sph(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const long num,
		    double semicol) {

  int n, m;
  assert(stride >0 && num >0);

  n = num * (stride + semicol) * d; // pre-allocate
  m = num * (stride + semicol); p_data_sph->max_col = m;

  p_data_sph->dim = d;  
  p_data_sph->str = stride;
  //  p_data_sph->size = num;


  p_data_sph->p_str  = _D2_CALLOC_INT(num);
  p_data_sph->p_str_cum  = _D2_CALLOC_LONG(num);
  p_data_sph->p_w    = _D2_MALLOC_SCALAR(m);

  if (d>0) 
    p_data_sph->p_supp = _D2_MALLOC_SCALAR(n);
  else if (d==0) 
    p_data_sph->dist_mat = _D2_MALLOC_SCALAR( stride * stride );


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

int d2_free(mph *p_data) {
  int i;
  for (i=0; i<p_data->s_ph; ++i) {
    if (p_data->ph[i].col > 0) d2_free_sph(p_data->ph + i);
  }
  free(p_data->ph);
  if (!p_data->label) _D2_FREE(p_data->label);
  return 0;
}

int d2_allocate_work(mph *p_data, var_mph *var_work) {
  int i;
  long size = p_data->size;
  int num_of_labels = p_data->num_of_labels;
  trieq *p_tr = &var_work->tr;
  var_work->s_ph = p_data->s_ph;

  var_work->g_var = (var_sph *) malloc(p_data->s_ph * sizeof(var_sph));
  if (d2_alg_type == D2_CENTROID_BADMM) 
    var_work->l_var_sphBregman = (var_sphBregman *) malloc(p_data->s_ph * sizeof(var_sphBregman));

  for (i=0; i<p_data->s_ph; ++i) {
    int str = p_data->ph[i].str;
    int col = p_data->ph[i].col;

    var_work->g_var[i].C = NULL;
    var_work->g_var[i].X = NULL;
    var_work->g_var[i].L = NULL;

    if (d2_alg_type == D2_CENTROID_BADMM || d2_alg_type == D2_CENTROID_ADMM) {
      var_work->g_var[i].C = _D2_MALLOC_SCALAR(str * col);
    }
    if (d2_alg_type == D2_CENTROID_BADMM) {
      d2_allocate_work_sphBregman(p_data->ph +i, p_data->size, 
				  var_work->l_var_sphBregman+i);
    }
    if (d2_alg_type == D2_CENTROID_ADMM) {
      var_work->g_var[i].X = _D2_MALLOC_SCALAR(str * col);
    }
    if (d2_alg_type == D2_CENTROID_GRADDEC || d2_alg_type == D2_CENTROID_ADMM) {
      var_work->g_var[i].X = _D2_MALLOC_SCALAR(str * col);
      var_work->g_var[i].L = _D2_MALLOC_SCALAR(str * size);
    }
  }

  var_work->label_switch = (char *) malloc(size * sizeof(char)); 

  p_tr->l = _D2_MALLOC_SCALAR(size * num_of_labels);
  p_tr->u = _D2_MALLOC_SCALAR(size);
  p_tr->s = _D2_MALLOC_SCALAR(num_of_labels);
  p_tr->c = _D2_MALLOC_SCALAR(num_of_labels * num_of_labels);
  p_tr->r = (char *) calloc(size, sizeof(char));

  return 0;
}

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
  free(var_work->g_var);
  if (d2_alg_type == D2_CENTROID_BADMM) free(var_work->l_var_sphBregman);
  free(var_work->label_switch);

  _D2_FREE(p_tr->l);
  _D2_FREE(p_tr->u);
  _D2_FREE(p_tr->s);
  _D2_FREE(p_tr->c);
  free(p_tr->r);
  return 0;
}

/** Compute the distance from each point to the all centroids.
    This part can be parallelized.
 */
long d2_labeling_prep(__IN_OUT__ mph *p_data,
		      mph *centroids,
		      var_mph * var_work,
		      int selected_phase) {
  long i, count = 0, dist_count = 0;
  int **p_str;
  double **p_supp, **p_w;
  long **p_str_cum;
  const long size = p_data->size;
  const int num_of_labels = centroids->size;
  int *label = p_data->label;
  trieq *p_tr = &var_work->tr;

  nclock_start();

  {
  int n;
  int s_ph = p_data->s_ph;

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
  }

  /* step 1 */
  for (i=0; i<num_of_labels; ++i) p_tr->s[i] = DBL_MAX;
#pragma omp parallel for reduction(+:dist_count)
  for (i=0; i<num_of_labels; ++i) {
    int n;
    long j;
    p_tr->c[i*num_of_labels + i] = 0; // d(c_i, c_i)

    for (j=i+1; j<num_of_labels; ++j) {
      double d = 0.0, val;
      for (n=0; n<p_data->s_ph; ++n) 
	if (selected_phase < 0 || n == selected_phase) {
	  int str = centroids->ph[n].str;
	  int dim = p_data->ph[n].dim;
	  assert(dim == centroids->ph[n].dim);
	  if (dim > 0) {
	  val = d2_match_by_coordinates(dim, 
					str, centroids->ph[n].p_supp + i*str*dim, centroids->ph[n].p_w + i*str, 
					str, centroids->ph[n].p_supp + j*str*dim, centroids->ph[n].p_w + j*str, 
					NULL, // x and lambda are implemented later
					NULL);
	  } else if (dim == 0) {
	    val = d2_match_by_distmat(str, str, centroids->ph[n].dist_mat, centroids->ph[n].p_w + i*str, centroids->ph[n].p_w + j*str, NULL, NULL);
	  }
	  d += val;
	}
      d = sqrt(d); dist_count +=1;

      p_tr->c[i*num_of_labels + j] = d; 
      p_tr->c[i + j*num_of_labels] = d;

      if (p_tr->s[i] > d) p_tr->s[i] = d;
      if (p_tr->s[j] > d) p_tr->s[j] = d;
    }    
  }


  for (i=0; i<size; ++i) 
    if (d2_alg_type == 0) {
      var_work->label_switch[i] = 0;
    }

#pragma omp parallel for reduction(+:dist_count,count)
  for (i=0; i<size; ++i) {
  /* step 2 */
  if (p_tr->u[i] > p_tr->s[label[i]]) {
    int init_label = label[i];
    int jj = init_label>=0? init_label: 0;
    int n;
    long j;
    SCALAR min_distance;
    SCALAR *U = p_tr->u + i;
    SCALAR *L = p_tr->l + i*num_of_labels;

    /* step 3 */
    for (j=0; j<num_of_labels; ++j) 
      if ((j != jj || init_label < 0) && *U > L[j] && *U > p_tr->c[j*num_of_labels + jj] / 2.) {

	/* 3a. */
	if (p_tr->r[i] == 1) {
	  /* compute distance */
	  double d = 0.0, val;
	  for (n=0; n<p_data->s_ph; ++n) 
	    if (selected_phase < 0 || n == selected_phase) {
	      int str = centroids->ph[n].str;
	      int dim = p_data->ph[n].dim;
	      assert(dim == centroids->ph[n].dim);
	      if (dim > 0) {
	      val = d2_match_by_coordinates(dim, 
					    p_str[n][i], p_supp[n] + p_str_cum[n][i]*dim, p_w[n] + p_str_cum[n][i],
					    str, centroids->ph[n].p_supp + jj*str*dim, centroids->ph[n].p_w + jj*str, 
					    NULL, // x and lambda are implemented later
					    NULL);
	      } else if (dim == 0) {
		val = d2_match_by_distmat(p_str[n][i], str, centroids->ph[n].dist_mat, p_w[n] + p_str_cum[n][i], centroids->ph[n].p_w + jj*str, NULL, NULL);
	      }

	      d += val;
	    }
	  d = sqrt(d); dist_count +=1;
	  L[jj] = d;
	  *U = d;
	  min_distance = d;
	  p_tr->r[i] = 0;
	} else {
	  min_distance = *U;
	}

	/* 3b. */
	if ((min_distance > L[j] || min_distance > p_tr->c[j*num_of_labels + jj] / 2.) && j!=jj) {
	  /* compute distance */
	  double d = 0.0, val;
	  for (n=0; n<p_data->s_ph; ++n) 
	    if (selected_phase < 0 || n == selected_phase) {
	      int str = centroids->ph[n].str;
	      int dim = p_data->ph[n].dim;
	      assert(dim == centroids->ph[n].dim);
	      if (dim > 0) {
	      val = d2_match_by_coordinates(dim, 
					    p_str[n][i], p_supp[n] + p_str_cum[n][i]*dim, p_w[n] + p_str_cum[n][i],
					    str, centroids->ph[n].p_supp + j*str*dim, centroids->ph[n].p_w + j*str, 
					    NULL, // x and lambda are implemented later
					    NULL);
	      } else if (dim == 0) {
		val = d2_match_by_distmat(p_str[n][i], str, centroids->ph[n].dist_mat, p_w[n] + p_str_cum[n][i], centroids->ph[n].p_w + j*str, NULL, NULL);
	      }

	      d += val;
	    }
	  d = sqrt(d); dist_count +=1;
	  L[j] = d;
	  if (d < min_distance) {jj = j; min_distance = d; *U = d;}
	}
      }
    
    if (jj != init_label) {
      label[i] = jj;
      if (d2_alg_type == 0) {
	var_work->label_switch[i] = 1;
      }
      count += 1;
    }
  }
  }

  free(p_str); free(p_supp); free(p_w); free(p_str_cum);
  VPRINTF(("\n\t\t\t\t %ld objects change their labels\n\t\t\t\t %ld distance pairs computed\n\t\t\t\t seconds: %f\n", count, dist_count, nclock_end()));
  
  return count;
}


/** copy a to b */
int d2_copy(mph* a, mph *b) {
  int n;
  bool new_init_tag = false;
  b->s_ph = a->s_ph;
  b->size = a->size;
  if (!b->ph) {
    b->ph = (sph *) malloc(b->s_ph * sizeof(sph));
    new_init_tag = true;
  }
  for (n=0; n<a->s_ph; ++n) 
    // check whether n-th phase is allocated
    if (a->ph[n].col > 0) {
      // check whether b->ph[n] is allocated; if not, allocate first
      if (new_init_tag)  {
	d2_allocate_sph(b->ph + n, a->ph[n].dim, a->ph[n].str, a->size, 0.);
	b->ph[n].col = a->ph[n].col;
      }
      memcpy(b->ph[n].p_str, a->ph[n].p_str, a->size * sizeof(int));
      memcpy(b->ph[n].p_str_cum, a->ph[n].p_str_cum, a->size * sizeof(long));
      memcpy(b->ph[n].p_w, a->ph[n].p_w, a->ph[n].col * sizeof(SCALAR));
      if (a->ph[n].dim > 0) 
	memcpy(b->ph[n].p_supp, a->ph[n].p_supp, a->ph[n].col * a->ph[n].dim * sizeof(SCALAR));
      else if (a->ph[n].dim == 0)  
	memcpy(b->ph[n].dist_mat, a->ph[n].dist_mat, a->ph[n].str*a->ph[n].str);
    } else {
      b->ph[n].col = 0;
    }
  
  return 0;
}

/** Compute the distance from each point to the all centroids.
    This part can be parallelized.
 */

#define max(X,Y) (((X) > (Y)) ? (X) : (Y))
long d2_labeling_post(mph *p_data,
		      mph *c_old,
		      mph *c_new,
		      var_mph * var_work,
		      int selected_phase) {
  int i, num_of_labels = c_old->size;
  long j, size = p_data->size;
  SCALAR *d_changes = _D2_MALLOC_SCALAR(num_of_labels);
  int *label = p_data->label;

  for (i=0; i<num_of_labels; ++i) {
    double d = 0, val;
    int n;
    for (n=0; n<c_old->s_ph; ++n) 
      if (selected_phase < 0 || n == selected_phase) {
	int str = c_old->ph[n].str;
	int dim = c_old->ph[n].dim;
	assert(dim == c_new->ph[n].dim);
	if (dim > 0) {
	val = d2_match_by_coordinates(dim, 
				      str, c_old->ph[n].p_supp + i*str*dim, c_old->ph[n].p_w + i*str, 
				      str, c_new->ph[n].p_supp + i*str*dim, c_new->ph[n].p_w + i*str, 
				      NULL, // x and lambda are implemented later
				      NULL);
	} else if (dim == 0) {
	  val = d2_match_by_distmat(str, str, c_old->ph[n].dist_mat, c_old->ph[n].p_w + i*str,c_new->ph[n].p_w + i*str, NULL, NULL);
	}

	d += val;
      }
    d = sqrt(d);    
    d_changes[i] = d;
  }

  for (j=0; j<size; ++j) {
    SCALAR * L = var_work->tr.l + j*num_of_labels;
    for (i=0; i<num_of_labels; ++i) L[i] = max(L[i] - d_changes[i], 0);
    var_work->tr.u[j] += d_changes[label[j]];
    var_work->tr.r[j] = 1;
  }

  _D2_FREE(d_changes);
  return 0;
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
  mph the_centroids_copy = (mph) {0, 0, NULL, 0, NULL};

  assert(num_of_clusters>0 && max_iter > 0);

  // label all objects as invalid numbers
  p_data->num_of_labels = num_of_clusters;
  for (i=0; i<size; ++i) label[i] = -1; // rand() % num_of_clusters;

  // initialize centroids from random
  centroids->s_ph = s_ph;
  centroids->size = num_of_clusters;
  centroids->ph = (sph *) malloc(s_ph * sizeof(sph));
  for (i=0; i<s_ph; ++i) 
    if (selected_phase < 0 || i == selected_phase) {
      d2_centroid_rands(p_data, i, centroids->ph + i);
    } else {
      centroids->ph[i].col = 0;
    }
  //d2_write(NULL, centroids); getchar();

  // allocate initialize auxiliary variables
  d2_allocate_work(p_data, &var_work);
  for (i=0; i<size * num_of_clusters; ++i) var_work.tr.l[i] = 0;
  for (i=0; i<size; ++i) {var_work.tr.u[i] = DBL_MAX; var_work.tr.r[i] = 1; }

  // start centroid-based clustering here
  d2_solver_setup();  
  for (iter=0; iter<max_iter; ++iter) {
    VPRINTF(("Round %d ... \n", iter));
    VPRINTF(("\tRe-labeling all instances ... ")); VFLUSH;
    label_change_count = d2_labeling_prep(p_data, centroids, &var_work, selected_phase);


    /* termination criterion */
    if (label_change_count < 0.005 * size) {
      VPRINTF(("Terminate!\n"));
      break;
    }

    /* make copies of centroids */
    d2_copy(centroids, &the_centroids_copy);

    VPRINTF(("\tUpdate centroids ... \n"));
    /* update centroids */
    for (i=0; i<s_ph; ++i) 
      if (selected_phase < 0 || i == selected_phase) {
	VPRINTF(("\t phase %ld: \n", i));            
      
	if (d2_alg_type == D2_CENTROID_BADMM) 
	  d2_centroid_sphBregman(p_data, &var_work, i, centroids->ph + i, centroids->ph + i);
	if (d2_alg_type == D2_CENTROID_GRADDEC)
	  d2_centroid_sphGradDecent(p_data, &var_work, i, centroids->ph + i, centroids->ph + i);
	if (d2_alg_type == D2_CENTROID_ADMM)
	  d2_centroid_sphADMM(p_data, &var_work, i, centroids->ph + i, centroids->ph + i);
      }

    /* post updates */
    d2_labeling_post(p_data, &the_centroids_copy, centroids, &var_work, selected_phase);
  }
  //d2_solver_release();

  d2_free_work(&var_work);
  d2_free(&the_centroids_copy);
  return 0;
}


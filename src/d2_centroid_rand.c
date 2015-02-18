#include "d2_clustering.h"
#include "d2_math.h"
#include "d2_param.h"
#include <stdbool.h>
#include <float.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/* initialize with multivariate normal */
/*
int d2_centroid_randn(mph *p_data, int idx_ph, sph *c) {
  int i;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  SCALAR *means, *covs;
    
  // set stride
  for (i=0; i<num_of_labels; ++i) {
    c->p_str[i] = str;
    c->p_str_cum[i] = i*str;
  }
  
  // set weight
  for (i=0; i<num_of_labels*str; ++i) {
    c->p_w[i] = (double) 1./(double) str;
  }

  // set column
  c->col = str * num_of_labels;

  // set dist_mat
  if (dim == 0) {
    for (i=0; i<str * str; ++i) c->dist_mat[i] = data_ph->dist_mat[i];
  }
  
  // compute mean and cov
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
*/
void merge_symbolic(const int dim,
		    const int * m_supp, const SCALAR * m_w, const int m,
		    int * c_supp, SCALAR * c_w, const int n,
		    const int vocab_size, 
		    const SCALAR* dist_mat) {
  SCALAR * D, *w;
  int *supp;
  int i, j, k, mm = m;

  assert(m>n);
  D = _D2_MALLOC_SCALAR(m*m);
  supp = _D2_MALLOC_INT(dim*m);
  w = _D2_MALLOC_SCALAR(m);

  for (i=0; i<dim*m; ++i) supp[i] = m_supp[i];
  for (i=0; i<m; ++i) w[i] = m_w[i];

  _D2_FUNC(pdist_symbolic)(dim, m, m, supp, supp, D, vocab_size, dist_mat);
  while (mm > n) {
    SCALAR min = DBL_MAX;
    int min_idx = -1;
    for (i=0; i<m; ++i)
      for (j=i+1; j<m; ++j)
	if (w[i] > 0 && w[j] > 0 && min > D[i*m+j]) {
	  min = D[i*m + j];
	  min_idx = i*m + j;
	}
    i = min_idx / m;
    j = min_idx % m;

    //merge (i,j)
    for (k=0; k<dim; ++k)
      supp[i*dim + k] = w[i] > w[j]? supp[i*dim + k]:supp[j*dim + k];

    w[i] = w[i] + w[j];
    w[j] = 0;
    
    mm-- ;
  }

  for (i=0, j=0; i<m; ++i) 
    if (w[i] > 0) {
      for (k=0; k<dim; ++k) c_supp[j*dim+k] = supp[i*dim+k];
      c_w[j] = w[i];
      j++;
    }
  assert(j==n);

  _D2_FREE(D); _D2_FREE(supp); _D2_FREE(w);
}

void merge         (const int dim, 
		    const SCALAR * m_supp, const SCALAR * m_w, const int m, 
		    SCALAR * c_supp, SCALAR * c_w, const int n) {
  SCALAR * D, *supp, *w;
  int i, j, k, mm = m;

  assert(m>n);

  D = _D2_MALLOC_SCALAR(m*m);
  supp = _D2_MALLOC_SCALAR(dim*m);
  w = _D2_MALLOC_SCALAR(m);
  
  _D2_CBLAS_FUNC(copy)(dim*m, m_supp, 1, supp, 1);
  _D2_CBLAS_FUNC(copy)(m, m_w, 1, w, 1);

  _D2_FUNC(pdist2)(dim, m, m, supp, supp, D);

  while (mm > n) {
    SCALAR min = DBL_MAX; // a very large number
    int min_idx = -1;
    for (i=0; i<m; ++i) 
      for (j=i+1; j<m; ++j) 
	if (w[i] > 0 && w[j] > 0 && min > D[i*m + j]) {
	  min = D[i*m + j];
	  min_idx = i*m + j;
	}
    i = min_idx / m;
    j = min_idx % m;
    // printf("%d %d\n", i, j);
    
    // merge (i,j)
    for (k=0; k < dim; ++k) 
      supp[i*dim + k] = (supp[i*dim + k] * w[i] + supp[j*dim + k] * w[j]) / (w[i] + w[j]);

    _D2_FUNC(pdist2)(dim, 1, m, supp + i*dim, supp, D + i*m);
    for (k=0; k < m; ++k) D[k*m + i] = D[i*m + k];

    w[i] = w[i] + w[j];
    w[j] = 0;
    
    mm-- ;
  }


  for (i=0, j=0; i<m; ++i) 
    if (w[i] > 0) {
      _D2_CBLAS_FUNC(copy)(dim, supp + i*dim, 1, c_supp + j*dim, 1);
      c_w[j] = w[i];
      j++;
    }
  assert(j == n);
  _D2_FREE(D); _D2_FREE(supp); _D2_FREE(w);
}

/* initialize with random samples */
int d2_centroid_rands(mph *p_data, int idx_ph, sph *c) {
  size_t i, j, k, *array;
  sph *data_ph = p_data->ph + idx_ph;
  size_t num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  int vocab_size = data_ph->vocab_size;
  size_t size = p_data->size;
  int strxdim = str*dim;

  SCALAR *m_supp, *m_w;
  int *m_supp_sym, d;

  assert(c->str == str);

  // set stride
  for (i=0; i<num_of_labels; ++i) {
    c->p_str[i] = str;
    c->p_str_cum[i] = i*str;
  }

  // set column
  c->col = str * num_of_labels;

  // set vocab_size and dist_mat
  if (data_ph->metric_type == D2_HISTOGRAM || data_ph->metric_type == D2_N_GRAM) {
    for (i=0; i<vocab_size * vocab_size; ++i) c->dist_mat[i] = data_ph->dist_mat[i];
    c->vocab_size = data_ph->vocab_size;
  }
  
  // generate index array
  array = _D2_MALLOC_SIZE_T(size);
  for (i = 0; i < size; ++i) array[i] = i;
  shuffle(array, size);

  // set to zero
  for (i=0; i<c->col; ++i) c->p_w[i] = 0;
  if (c->metric_type == D2_EUCLIDEAN_L2) {
    for (i=0; i<c->col * c->dim; ++i) c->p_supp[i] = 0;
  }
  else if (c->metric_type == D2_N_GRAM) {
    for (i=0; i<c->col * c->dim; ++i) c->p_supp_sym[i] = 0;
  }

  

  i = 0; j = world_rank;

  while (i<size && j<num_of_labels) {
    int the_str;
    size_t the_str_cum;

    while (i<size && data_ph->p_str[array[i]] < str)  { ++i; }
    if (i == size) break;
    the_str = data_ph->p_str[array[i]];
    the_str_cum = data_ph->p_str_cum[array[i]];

    switch (data_ph->metric_type) {
    case D2_EUCLIDEAN_L2:
    case D2_HISTOGRAM:
      m_supp = data_ph->p_supp + the_str_cum*dim; 
      m_w = data_ph->p_w + the_str_cum;

      if (the_str == str) {
	_D2_CBLAS_FUNC(copy)(strxdim, m_supp, 1, c->p_supp + j*strxdim , 1);
	_D2_CBLAS_FUNC(copy)(str, m_w, 1, c->p_w + j*str , 1);
      } else {
	merge(dim, 
	      m_supp, m_w, the_str, 
	      c->p_supp + j*strxdim, c->p_w + j*str, str);
      }
      break;
    case D2_WORD_EMBED:
      m_supp_sym = data_ph->p_supp_sym + the_str_cum;
      m_w = data_ph->p_w + the_str_cum;
      if (the_str == str) {
	for (k=0; k<str; ++k) {
	  for (d=0; d<dim; ++d) 
	    c->p_supp[j*strxdim + k*dim + d] = data_ph->vocab_vec[m_supp_sym[k]*dim + d];
	  c->p_w[j*str + k] = m_w[k];
	}
      } else {
	SCALAR *supp = (SCALAR*) _D2_MALLOC_SCALAR(the_str*dim);
	for (k=0; k<the_str; ++k)
	  for (d=0; d<dim; ++d)
	    supp[k*dim + d] = data_ph->vocab_vec[m_supp_sym[k]*dim + d];
	merge(dim,
	      supp, m_w, the_str,
	      c->p_supp + j*strxdim, c->p_w + j*str, str);
	_D2_FREE(supp);
      }
      break;
    case D2_N_GRAM: 
      m_supp_sym = data_ph->p_supp_sym + the_str_cum*dim;
      m_w = data_ph->p_w + the_str_cum;
      if (the_str == str) {
	for (k=0; k<strxdim; ++k) c->p_supp_sym[j*strxdim + k] = m_supp_sym[k];
	for (k=0; k<str; ++k) c->p_w[j*str + k] = m_w[k];
      } else {
	merge_symbolic(dim, 
		       m_supp_sym, m_w, the_str, 
		       c->p_supp_sym + j*strxdim, c->p_w + j*str, str,
		       c->vocab_size, c->dist_mat);
      }
      break;
    default:
      fprintf(stderr, "Unknown type of data!");
      assert(false);
    }
    ++i; j+= nprocs;
  }
  assert(j>=num_of_labels);

  free(array);
  return 0;
}

#include "d2_clustering.h"
#include "d2_math.h"

/* initialize with multivariate normal */
int d2_centroid_randn(mph *p_data, int idx_ph, sph *c) {
  int i;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  SCALAR *means, *covs;
  
  // initialization 
  d2_allocate_sph(c, dim, str, num_of_labels, 0.);
  
  // set stride
  for (i=0; i<num_of_labels; ++i) {
    c->p_str[i] = str;
  }
  
  // set weight
  for (i=0; i<num_of_labels*str; ++i) {
    c->p_w[i] = (double) 1./(double) str;
  }

  // set column
  c->col = str * num_of_labels;
  
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


void merge(const int dim, const SCALAR * m_supp, const SCALAR * m_w, const int m, SCALAR * c_supp, SCALAR * c_w, const int n) {
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
    SCALAR min = 1E10; // a very large number
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
  int i, j, k, *array;
  sph *data_ph = p_data->ph + idx_ph;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int str = data_ph->str;
  int size = p_data->size;
  int strxdim = str*dim;
  SCALAR *m_supp, *m_w;

  // initialization 
  d2_allocate_sph(c, dim, str, num_of_labels, 0.);

  // set stride
  for (i=0; i<num_of_labels; ++i) {
    c->p_str[i] = str;
  }

  // set column
  c->col = str * num_of_labels;
  
  // generate index array
  array = (int *) malloc(size * sizeof(int));
  shuffle(array, size);

  i = 0; j = 0;
  m_supp = data_ph->p_supp; m_w = data_ph->p_w;
  while (i<size && j<num_of_labels) {
    while (i<size && data_ph->p_str[i] < str) 
      { m_supp = m_supp + data_ph->p_str[i]*dim; m_w = m_w + data_ph->p_str[i];  ++i; }
    if (i == size) break;

    if (data_ph->p_str[i] == str) {
      _D2_CBLAS_FUNC(copy)(strxdim, m_supp, 1, c->p_supp + j*strxdim , 1);
      _D2_CBLAS_FUNC(copy)(str, m_w, 1, c->p_w + j*str , 1);
    } else {
      merge(dim, m_supp, m_w, data_ph->p_str[i], c->p_supp + j*strxdim, c->p_w + j*str, str);
    }

    m_supp = m_supp + data_ph->p_str[i]*dim; m_w = m_w + data_ph->p_str[i];  ++i; ++j;
  }

  assert(j == num_of_labels);

  free(array);
  return 0;
}

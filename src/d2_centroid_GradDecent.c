#include "d2_clustering.h"
#include "d2_math.h"
#include "stdio.h"
#include "d2_param.h"
#include "d2_solver.h"
#include <float.h>
#include <assert.h>
#include "d2_centroid_util.h"

static GRADDEC_options graddec_default_options = {.maxIters = 5, .stepSize = 0.5};

GRADDEC_options *p_graddec_options = &graddec_default_options;

extern double time_budget;


int d2_centroid_sphGradDecent(mph *p_data,
			      var_mph * var_work,
			      int idx_ph,
			      sph *c0,
			      __OUT__ sph *c) {
  

  size_t i;
  int iter, nIter = p_graddec_options->maxIters;
  double fval0, fval = DBL_MAX;
  int *label = p_data->label;
  size_t num_of_labels = p_data->num_of_labels;
  int str, strxdim;
  size_t size = p_data->size;
  sph *data_ph = p_data->ph + idx_ph;
  int dim = data_ph->dim;
  int *p_str = data_ph->p_str;
  SCALAR *p_supp = data_ph->p_supp;
  int *p_supp_sym = data_ph->p_supp_sym;
  SCALAR *p_w = data_ph->p_w;
  size_t *p_str_cum = data_ph->p_str_cum;
  SCALAR *X = var_work->g_var[idx_ph].X;
  SCALAR *L = var_work->g_var[idx_ph].L;
  SCALAR *C = var_work->g_var[idx_ph].C;
  size_t *label_count;
  SCALAR *p_grad, *Zr=NULL, tmp;
  double step_size = p_graddec_options->stepSize;
  double startTime;

  assert(dim > 0); // current only support the D2 format

  if (!c0) {
    d2_centroid_rands(p_data, idx_ph, c);// For profile purpose only
  } else {
    *c = *c0; // warm start (for clustering purpose)
  }  
  str = c->str; assert(str > 0);
  strxdim = str * dim;

  /** Calculate labels counts:
      it could be possible that some clusters might not 
      have any instances
   */
  label_count = _D2_MALLOC_SIZE_T(num_of_labels);    
  for (i=0; i<num_of_labels; ++i) label_count[i] = 0;
  for (i=0; i<size; ++i) ++label_count[label[i]];
  for (i=0; i<num_of_labels; ++i) assert(label_count[i] != 0);
  

  // allocate mem for gradient of p_w
  p_grad = _D2_MALLOC_SCALAR(num_of_labels * str);

  // conditonal allocation
  if (data_ph->metric_type == D2_N_GRAM) 
    Zr = _D2_MALLOC_SCALAR(str * num_of_labels * (strxdim * data_ph->vocab_size + 1));

  /* compute exact distances */
  calculate_distmat(data_ph, label, size, c, C);
  startTime = getRealTime();  
  for (iter = 0; iter <= nIter; ++iter) {

    fval0 = fval;
    fval = 0;

    for (i=0; i<size; ++i) {
      fval += d2_match_by_distmat(str, p_str[i],
				  C + str*p_str_cum[i],
				  c->p_w + label[i]*str, p_w + p_str_cum[i],
				  X + p_str_cum[i]*str,
				  L + i*str,
				  i); // known bug as work for only one phase      
    }
    fval /= size;
  
    printf("\t%d\t%f\t%f\n", iter, fval, getRealTime() - startTime);    
    if (fabs(fval - fval0) < 1E-4 * fval || getRealTime() - startTime > time_budget) break;

    /* compute p_grad */
    _D2_FUNC(ccenter)(str, size, L, NULL);
    for (i=0; i<num_of_labels*str; ++i) p_grad[i] = 0; //reset
    for (i=0; i < size; ++i)
      _D2_CBLAS_FUNC(axpy)(str, 1./label_count[label[i]], L + i*str, 1, p_grad + label[i]*str, 1);
    for (i=0, tmp = 0.; i < num_of_labels*str; ++i) {
      p_grad[i] = p_grad[i] * c->p_w[i] * (1-c->p_w[i]);
      tmp += p_grad[i]*p_grad[i];
    }

    /* update c->p_supp */
    switch (data_ph->metric_type) {
    case D2_EUCLIDEAN_L2 :
      for (i=0; i<strxdim*num_of_labels; ++i) c->p_supp[i] = 0; // reset c->p_supp
      for (i=0; i < size;  ++i) {
	_D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasTrans, 
			     dim, str, p_str[i], 
			     1, 
			     p_supp + dim*p_str_cum[i], dim, 
			     X + str*p_str_cum[i], str, 
			     1, 
			     c->p_supp + label[i]*strxdim, dim);
      }
      for (i=0; i<num_of_labels; ++i) {
	_D2_CBLAS_FUNC(scal)(str*dim, 1./label_count[i], c->p_supp + i*strxdim, 1);
	_D2_FUNC(irms)(dim, str, c->p_supp + i*strxdim, c->p_w + i*str);
      }
    
      /* compute exact distances */
      calculate_distmat(data_ph, label, size, c, C);

      break;
    case D2_WORD_EMBED :
      assert(num_of_labels < size);
      for (i=0; i<strxdim*num_of_labels; ++i) c->p_supp[i] = 0.f;
      for (i=0; i<size; ++i) {
	int *m_supp_sym = p_supp_sym + p_str_cum[i];
	SCALAR *c_supp = c->p_supp + label[i]*strxdim;
	SCALAR *Xm = X + str*p_str_cum[i];
	int j, k, d;

	for (k=0; k<str; ++k)
	  for (j=0; j<p_str[i]; ++j) 
	    if (Xm[j*str + k] > 1.E-9) {
	      for (d=0; d<dim; ++d)
		c_supp[k*dim + d] += Xm[j*str + k] * data_ph->vocab_vec[m_supp_sym[j]*dim + d];
	    }
      }
      for (i=0; i<num_of_labels; ++i) {
	_D2_CBLAS_FUNC(scal)(str*dim, 1./label_count[i], c->p_supp + i*strxdim, 1);
	_D2_FUNC(irms)(dim, str, c->p_supp + i*strxdim, c->p_w + i*str);
      }

      /* compute exact distances */
      calculate_distmat(data_ph, label, size, c, C);

      break;
    case D2_N_GRAM :
      if (iter == nIter) {
      assert(num_of_labels * (strxdim * data_ph->vocab_size + 1) <= size);

      for (i=0; i<num_of_labels*strxdim*data_ph->vocab_size; ++i) Zr[i] = 0; //reset Zr to temporarily storage
      for (i=0; i<size; ++i) {
	accumulate_symbolic(dim, str, p_str[i], 
			    p_supp_sym + dim*p_str_cum[i], 
			    X + str*p_str_cum[i], 
			    Zr + label[i]*strxdim*data_ph->vocab_size,
			    data_ph->vocab_size);	  
      }

      for (i=0; i<num_of_labels; ++i) {
	minimize_symbolic(dim, str, c->p_supp_sym + i*strxdim, Zr + i*strxdim*data_ph->vocab_size, data_ph->vocab_size, data_ph->dist_mat, Zr + num_of_labels*strxdim*data_ph->vocab_size);
      }

      /* compute exact distances */
      calculate_distmat(data_ph, label, size, c, C);

      }
      break;    
    }

    /* update c->p_w */
    if (step_size / sqrt(tmp) < 5.) {
      tmp = step_size / sqrt(tmp);
    } else {
      tmp = 5.;
    }

    //    printf("%lf\n", tmp);

    for (i=0; i<num_of_labels*str; ++i) c->p_w[i] *= exp(- tmp * p_grad[i]);
    _D2_FUNC(cnorm)(str, num_of_labels, c->p_w, NULL);

  }

  _D2_FREE(label_count);
  _D2_FREE(p_grad); 
  if (Zr) _D2_FREE(Zr);
  return 0;
}

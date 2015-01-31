#include "d2_clustering.h"
#include "d2_math.h"
#include "stdio.h"
#include "d2_param.h"
#include "d2_centroid_util.h"
#include <float.h>
#include <assert.h>

// #define __USE_MPI__ /* macro for the use of MPI */


/* choose options */

static BADMM_options badmm_clu_options = {.maxIters = 100, .rhoCoeff = 5.f, .updatePerLoops = 10};
static BADMM_options badmm_cen_options = {.maxIters = 2000, .rhoCoeff = 1.f, .updatePerLoops = 10};


BADMM_options *p_badmm_options = &badmm_clu_options;
 
int d2_allocate_work_sphBregman(sph *ph, size_t size, var_sphBregman * var_phwork) {
  size_t i; int j;
  SCALAR tmp;
  var_phwork->X = _D2_MALLOC_SCALAR (ph->str * ph->col);
  var_phwork->Z = _D2_MALLOC_SCALAR (ph->str * ph->col);
  var_phwork->Xc= _D2_MALLOC_SCALAR (ph->col);
  var_phwork->Zr= _D2_MALLOC_SCALAR (ph->str * size);
  var_phwork->Y = _D2_CALLOC_SCALAR (ph->str * ph->col); // initialized

  for (i=0; i<size; ++i) {
    SCALAR *p_scal = var_phwork->Z + ph->str*ph->p_str_cum[i];
    tmp = 1./(ph->str * ph->p_str[i]);
    for (j=0; j<ph->str*ph->p_str[i]; ++j, ++p_scal) *p_scal = tmp;
  }
  return 0;
}

int d2_free_work_sphBregman(var_sphBregman *var_phwork) {
  _D2_FREE(var_phwork->X);
  _D2_FREE(var_phwork->Z);
  _D2_FREE(var_phwork->Y);
  _D2_FREE(var_phwork->Xc);
  _D2_FREE(var_phwork->Zr);
  return 0;
}


/**
 * See matlab/centroid_sphBregman.m 
 * for a prototype implementation in Matlab. 
 */
int d2_centroid_sphBregman(mph *p_data, /* local data */
			   var_mph * var_work,
			   int idx_ph, /* index of phases */
			   sph *c0, /* initial gauss for the centroid */
			   __OUT__ sph *c
			   ) {
  sph *data_ph = p_data->ph + idx_ph;
  int *label = p_data->label;
  size_t num_of_labels = p_data->num_of_labels;
  char *label_switch = var_work->label_switch;
  int dim = data_ph->dim;
  size_t col = data_ph->col;
  int str, strxdim;
  size_t size = p_data->size;
  int *p_str = data_ph->p_str;
  SCALAR *p_supp = data_ph->p_supp;
  SCALAR *p_w = data_ph->p_w;
  size_t *p_str_cum = data_ph->p_str_cum;
  int *p_supp_sym = data_ph->p_supp_sym;
  SCALAR *C = var_work->g_var[idx_ph].C;
  SCALAR *X = var_work->l_var_sphBregman[idx_ph].X;
  SCALAR *Y = var_work->l_var_sphBregman[idx_ph].Y;
  SCALAR *Z = var_work->l_var_sphBregman[idx_ph].Z;
  SCALAR *Xc= var_work->l_var_sphBregman[idx_ph].Xc;
  SCALAR *Zr= var_work->l_var_sphBregman[idx_ph].Zr; 

  /**
   * MPI notes: vector needs synchronized __USE_MPI__ : 
   * @param(Zr, c->p_w)
   * @param(c->p_supp/c->p_supp_sym) for metric_type D2_EUCLIDEAN_L2 and D2_N_GRAM
   * @param(primres, dualres) 
   */
  


  size_t i; int j;
  int max_niter = p_badmm_options->maxIters, iter;
  SCALAR rho, obj, primres, dualres;
  SCALAR tmp, *Z0;
  size_t *label_count;

  /* Initialization */
  if (!c0) {
    // MPI note: to be done only on one node
    d2_centroid_rands(p_data, idx_ph, c);// For profile purpose only!
#ifdef __USE_MPI__
    /* initialize c from one node, and broadcast to other nodes */
#endif
  } else {
    *c = *c0; // warm start (for clustering purpose)
  }  
  str = c->str; assert(str > 0);
  strxdim = str * dim;

  /* compute C */  
  calculate_distmat(data_ph, label, size, c, C);

  /* rho is an important hyper-parameter */
  rho = p_badmm_options->rhoCoeff * _D2_CBLAS_FUNC(asum)(str*col, C, 1) / (str*col);

  if (dim > 0) { // for histogram, no re-init needed
  /* 
   * Indeed, we may only need to reinitialize for entries 
   * whose label are changed  
   */
    for (i=0; i<size; ++i) {
      if (label_switch[i] == 1) {
	SCALAR *p_scal = Z + str*p_str_cum[i];
	tmp = 1./(str * p_str[i]);
	for (j=0; j<str*p_str[i]; ++j, ++p_scal) *p_scal = tmp;
      }
    }
  }
  // allocate buffer of Z
  Z0 = _D2_MALLOC_SCALAR(str*col);

  /**
   *  Calculate labels counts:
   *  it could be possible that some clusters might not have any instances,
   *  thus this part needs improvement
   */
  label_count = _D2_CALLOC_SIZE_T(num_of_labels);    
  for (i=0; i<size; ++i) ++label_count[label[i]];
  for (i=0; i<num_of_labels; ++i) assert(label_count[i] != 0);

  // main loop
  printf("\titer\tobj\t\tprimres\t\tdualres\t\tseconds\n");
  printf("\t----------------------------------------------------------------\n");
  nclock_start();
  for (iter=0; iter <= max_niter; ++iter) {
    /*************************************************************************/
    // step 1: update X

    // X = Z.*exp(- (C + Y)/rho)
    _D2_CBLAS_FUNC(copy)(str*col, C, 1, X, 1);
    _D2_CBLAS_FUNC(axpy)(str*col, 1, Y, 1, X, 1);
    _D2_CBLAS_FUNC(scal)(str*col, -1./rho, X, 1);
    _D2_FUNC(exp)(str*col, X);
    _D2_FUNC(vmul)(str*col, X, Z, X);
    // normalize X
    _D2_FUNC(cnorm)(str, col, X, Xc); 
    _D2_FUNC(grms)(str, col, X, p_w);

    /*************************************************************************/
    // step 2: update Z
    // Z = X.*exp(Y/rho)
    _D2_CBLAS_FUNC(copy)(str*col, Z, 1, Z0, 1);
    _D2_CBLAS_FUNC(copy)(str*col, Y, 1, Z, 1);
    _D2_CBLAS_FUNC(scal)(str*col, 1./rho, Z, 1);
    _D2_FUNC(exp)(str*col, Z);
    _D2_FUNC(vmul)(str*col, Z, X, Z);
    
    // normalize Z
    for (i=0;i<size; ++i) {
      _D2_FUNC(rnorm)(str, p_str[i], Z + str*p_str_cum[i], Zr + str*i); 
      _D2_FUNC(gcms)(str, p_str[i], Z + str*p_str_cum[i], c->p_w + str*label[i]);
    }
   
    /*************************************************************************/
    // step 3: update Y
    _D2_CBLAS_FUNC(axpy)(str*col, rho, X, 1, Y, 1);
    _D2_CBLAS_FUNC(axpy)(str*col, -rho, Z, 1, Y, 1);


    /*************************************************************************/
    // step 4: update c->p_w
    for (i=0; i<str*num_of_labels; ++i) c->p_w[i] = 0; //reset c->p_w locally
    _D2_FUNC(cnorm)(str, size, Zr, NULL);


    // ADD vec(&Zr[str*i], str) TO vec(&c->p_w[label[i]*str], str)
    for (i=0; i<size; ++i) 
      _D2_CBLAS_FUNC(axpy)(str, 1, Zr + str*i, 1, c->p_w +label[i]*str, 1);
#ifdef __USE_MPI__
    /* ALLREDUCE by SUM operator: vec(c->p_w, num_of_labels*str) */
#endif
    _D2_FUNC(cnorm)(str, num_of_labels, c->p_w, NULL);
    

    // step 5: update c->p_supp (optional)
    if (iter % p_badmm_options->updatePerLoops == 0) {
      switch (data_ph->metric_type) {
      case D2_EUCLIDEAN_L2 :
	for (i=0; i<strxdim*num_of_labels; ++i) c->p_supp[i] = 0; // reset c->p_supp
	for (i=0; i<str*num_of_labels; ++i) Zr[i] = 0; //reset Zr to temporarily storage

	for (i=0;i < size;  ++i) {
	  /* ADD mat(&p_supp[p_str_cum[i]*dim], dim, p_str[i]) * 
	         mat(&X[p_str_cum[i]*str], str, p_str[i]).transpose 
	     TO mat(&c->p_supp[label[i]*strxdim], dim, str)
	   */
	  _D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasTrans, 
			       dim, str, p_str[i], 1, p_supp + dim*p_str_cum[i], dim, X + str*p_str_cum[i], str, 1, 
			       c->p_supp + label[i]*strxdim, dim);
	  /* ADD row_of_sums of mat(&X[p_str_cum[i]*str], str, p_str[i]) 
	     To vec(&Zr[label[i]*str], str) */
	  _D2_FUNC(rsum2)(str, p_str[i], X + str*p_str_cum[i], Zr + label[i]*str);
	}
#ifdef __USE_MPI__
	/* ALLREDUCE by SUM operator: vec(c->p_supp, num_of_labels*str*dim) */
	/* ALLREDUCE by SUM operator: vec(Zr, num_of_labels*str) */
#endif
	for (i=0; i<num_of_labels; ++i) {
	  _D2_FUNC(irms)(dim, str, c->p_supp + i*strxdim, Zr + i*str);
	}


	// re-calculate C
	calculate_distmat(data_ph, label, size, c, C);
	break;

      case D2_HISTOGRAM :
	break;

      case D2_N_GRAM :
	if (iter > 0) {
	assert(num_of_labels * (strxdim * data_ph->vocab_size + 1) <= size);

	for (i=0; i<num_of_labels*strxdim*data_ph->vocab_size; ++i) Zr[i] = 0; //reset Zr to temporarily storage
	for (i=0; i<size; ++i) {
	  accumulate_symbolic(dim, str, p_str[i], 
			      p_supp_sym + dim*p_str_cum[i], 
			      X + str*p_str_cum[i], 
			      Zr + label[i]*strxdim*data_ph->vocab_size,
			      data_ph->vocab_size);	  
	}
#ifdef __USE_MPI__
	/* ALLREDUCE by SUM operator: vec(Zr, num_of_labels*str*dim*data_ph->vocab_size) */
#endif
	for (i=0; i<num_of_labels; ++i) {
	  minimize_symbolic(dim, str, c->p_supp_sym + i*strxdim, Zr + i*strxdim*data_ph->vocab_size, data_ph->vocab_size, data_ph->dist_mat, Zr + num_of_labels*strxdim*data_ph->vocab_size);
	}

	// re-calculate C
	calculate_distmat(data_ph, label, size, c, C);
	}
	break;
      }
    }    

    /*************************************************************************/
    // step 6: check residuals
    obj = _D2_CBLAS_FUNC(dot)(str*col, C, 1, X, 1) / size;
    if (iter%100 == 0 || (iter < 100 && iter%20 == 0) ) {
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, X, 1);
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, Z0,1);
      primres = _D2_CBLAS_FUNC(asum)(str*col, X, 1) / size;
      dualres = _D2_CBLAS_FUNC(asum)(str*col,Z0, 1) / size;
#ifdef __USE_MPI__
      /* ALLREDUCE by MEAN operator: primres, dualres */
#endif
      printf("\t%d\t%f\t%f\t%f\t%f\n", iter, obj, primres, dualres, nclock_end());
    }
  }

  _D2_FREE(Z0);
  _D2_FREE(label_count);

  return 0;
}

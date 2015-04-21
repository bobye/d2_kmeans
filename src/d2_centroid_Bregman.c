#include "d2_clustering.h"
#include "d2_math.h"
#include "stdio.h"
#include "d2_param.h"
#include "d2_centroid_util.h"
#include <float.h>
#include <assert.h>

#ifdef __USE_MPI__
#include <mpi.h>
#endif



/* choose options */

static BADMM_options badmm_clu_options = {.maxIters = 100, .rhoCoeff = 2.f, .updatePerLoops = 10};
static BADMM_options badmm_cen_options = {.maxIters = 2000, .rhoCoeff = 1.f, .updatePerLoops = 10};

#define ROUNDOFF (1E-9)

BADMM_options *p_badmm_options = &badmm_clu_options;

extern double time_budget;
 
int d2_allocate_work_sphBregman(sph *ph, size_t size, var_sphBregman * var_phwork) {
  assert(ph->str > 0 && ph->col > 0 && size > 0);
  var_phwork->X = _D2_MALLOC_SCALAR (ph->str * ph->col); assert(var_phwork->X);
  var_phwork->Z = _D2_MALLOC_SCALAR (ph->str * ph->col); assert(var_phwork->Z);
  var_phwork->Xc= _D2_MALLOC_SCALAR (ph->col);           assert(var_phwork->Xc);
  var_phwork->Zr= _D2_MALLOC_SCALAR (ph->str * size);    assert(var_phwork->Zr);
  var_phwork->Y = _D2_MALLOC_SCALAR (ph->str * ph->col); assert(var_phwork->Y); // initialized

  return 0;
}

int d2_free_work_sphBregman(var_sphBregman *var_phwork) {
  if (var_phwork->X) _D2_FREE(var_phwork->X);
  if (var_phwork->Z) _D2_FREE(var_phwork->Z);
  if (var_phwork->Y) _D2_FREE(var_phwork->Y);
  if (var_phwork->Xc)_D2_FREE(var_phwork->Xc);
  if (var_phwork->Zr)_D2_FREE(var_phwork->Zr);
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
  SCALAR *Zr2 = Zr; char hasZr2 = 0;

  /**
   * MPI notes: vector needs synchronized __USE_MPI__ : 
   * @param(Zr, c->p_w)
   * @param(c->p_supp/c->p_supp_sym) for metric_type D2_EUCLIDEAN_L2 and D2_N_GRAM
   * @param(primres, dualres) 
   */
  


  size_t i; int j;
  int max_niter = p_badmm_options->maxIters, iter;
  SCALAR rho, obj, primres, dualres;
  SCALAR *Z0;
  size_t *label_count;

  /* Initialization */
  if (!c0) {
    // MPI note: to be done only on one node
    d2_centroid_rands(p_data, idx_ph, c);// For profile purpose only!
    broadcast_centroids(p_data, idx_ph);
  } else {
    *c = *c0; // warm start (for clustering purpose)
  }  
  str = c->str; assert(str > 0);
  strxdim = str * dim;

  /* compute C */  
  calculate_distmat(data_ph, label, size, c, C);

  /* rho is an important hyper-parameter */
  rho = p_badmm_options->rhoCoeff * _D2_CBLAS_FUNC(asum)(str*col, C, 1) / (str*col);
  for (i=0; i<str*col; ++i) C[i] /= rho; // normalize C and Y

  if (dim > 0) { // for histogram, no re-init needed
  /* 
   * Indeed, we may only need to reinitialize for entries 
   * whose label are changed  
   */
    for (i=0; i<size; ++i) {
      if (label_switch[i] == 1) {
	SCALAR *p_scal = Z + str*p_str_cum[i];
	SCALAR *data_w_scal  = p_w + p_str_cum[i];
	SCALAR *c_w_scal = c->p_w + label[i]*str;
	for (j=0; j<str*p_str[i]; ++j, ++p_scal) 
	  *p_scal =  data_w_scal[j/str] * c_w_scal[j%str];
      }
    }
  }
  // allocate buffer of Z
  Z0 = _D2_MALLOC_SCALAR(str*col);
  for (i=0; i<str*col; ++i) Y[i] = 0; // set Y to zero
  if (str * num_of_labels * (strxdim * data_ph->vocab_size + 1) > size*str - c->col && data_ph->metric_type == D2_N_GRAM) {
    Zr2 = _D2_MALLOC_SCALAR(str * num_of_labels * (strxdim * data_ph->vocab_size + 1));    
    hasZr2 = 1;
  }
  /**
   *  Calculate labels counts:
   *  it could be possible that some clusters might not have any instances,
   *  thus this part needs improvement
   */
  label_count = _D2_CALLOC_SIZE_T(num_of_labels);    
  for (i=0; i<size; ++i) ++label_count[label[i]];
#ifdef __USE_MPI__
  assert(sizeof(size_t) == sizeof(unsigned long long));
  MPI_Allreduce(MPI_IN_PLACE, label_count, num_of_labels, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif  
  VPRINTF("\tlabel counts:"); for (i=0; i<num_of_labels; ++i) {assert(label_count[i] != 0); VPRINTF("%d ", label_count[i]);} VPRINTF("\n");

  // main loop
  VPRINTF("\titer\tobj\t\tprimres\t\tdualres\t\tseconds\n");
  VPRINTF("\t----------------------------------------------------------------\n");
  nclock_start();
  for (iter=0; iter <= max_niter; ++iter) {
    /*************************************************************************/
    // step 1: update X

    // X = Z.*exp(- (C + Y)/rho)
    for (i=0; i<str*col; ++i) {
      X[i] = Z[i] * exp (- (C[i] + Y[i])) + ROUNDOFF;
    }      
    _D2_FUNC(cnorm)(str, col, X, Xc); 
    _D2_FUNC(grms)(str, col, X, p_w);

    /*************************************************************************/
    // step 2: update Z
    // Z = X.*exp(Y/rho)
    for (i=0; i<str*col; ++i) {
      Z0[i] = Z[i];
      Z[i]  = X[i] * exp(Y[i]) + ROUNDOFF;
    }
    for (i=0;i<size; ++i) {
      _D2_FUNC(rnorm)(str, p_str[i], Z + str*p_str_cum[i], Zr + str*i); 
      _D2_FUNC(gcms)(str, p_str[i], Z + str*p_str_cum[i], c->p_w + str*label[i]);
    }
   
    /*************************************************************************/
    // step 3: update Y
    for (i=0; i<str*col; ++i) Y[i] += X[i] - Z[i];

    /*************************************************************************/
    // step 4: update c->p_w
    for (i=0; i<str*num_of_labels; ++i) c->p_w[i] = 0; //reset c->p_w locally
    _D2_FUNC(cnorm)(str, size, Zr, Xc);


    // ADD vec(&Zr[str*i], str) TO vec(&c->p_w[label[i]*str], str)
    for (i=0; i<size; ++i) 
      _D2_CBLAS_FUNC(axpy)(str, 1, Zr + str*i, 1, c->p_w +label[i]*str, 1);
#ifdef __USE_MPI__
    /* ALLREDUCE by SUM operator: vec(c->p_w, c->col) */
    MPI_Allreduce(MPI_IN_PLACE, c->p_w, c->col, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    _D2_FUNC(cnorm)(str, num_of_labels, c->p_w, Xc);
    

    // step 5: update c->p_supp (optional)
    if (iter % p_badmm_options->updatePerLoops == 0) {
      switch (data_ph->metric_type) {
      case D2_EUCLIDEAN_L2 :
	//	assert(num_of_labels < size);
	for (i=0; i<strxdim*num_of_labels; ++i) c->p_supp[i] = 0.f; // reset c->p_supp
	for (i=0; i<c->col; ++i) Zr[i] = 0.f; //reset Zr to temporarily storage

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
	/* ALLREDUCE by SUM operator: vec(c->p_supp, c->col*dim) */
	MPI_Allreduce(MPI_IN_PLACE, c->p_supp, c->col * dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	/* ALLREDUCE by SUM operator: vec(Zr, c->col) */
	MPI_Allreduce(MPI_IN_PLACE, Zr, c->col, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	for (i=0; i<num_of_labels; ++i) {
	  _D2_FUNC(irms)(dim, str, c->p_supp + i*strxdim, Zr + i*str);
	}


	// re-calculate C
	calculate_distmat(data_ph, label, size, c, C);
	/* rho is an important hyper-parameter */
	for (i=0; i<str*col; ++i) C[i] /= rho; // normalize C and Y
	break;
      case D2_WORD_EMBED :
	assert(num_of_labels < size);
	for (i=0; i<strxdim*num_of_labels; ++i) c->p_supp[i] = 0.f;
	for (i=0; i<c->col; ++i) Zr[i] = 0.f; //reset Zr to temporarily storage
	for (i=0; i<size; ++i) {
	  int *m_supp_sym = p_supp_sym + p_str_cum[i];
	  SCALAR *Xm = X + str*p_str_cum[i];
	  SCALAR *Zrm = Zr + label[i]*str;
	  int j, k, d;
	  for (j=0; j<p_str[i]; ++j) {
	    SCALAR *c_supp = c->p_supp + label[i]*strxdim;
	    SCALAR *vocab_vec=&data_ph->vocab_vec[m_supp_sym[j]*dim];
	    for (k=0; k<str; ++k, ++Xm) {
	      for (d=0; d<dim; ++d, ++c_supp)		
		*c_supp += *Xm *  vocab_vec[d];
	      Zrm[k] += *Xm;
	    }
	  }
	}
#ifdef __USE_MPI__
	/* ALLREDUCE by SUM operator: vec(c->p_supp, c->col*dim) */
	MPI_Allreduce(MPI_IN_PLACE, c->p_supp, c->col * dim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	/* ALLREDUCE by SUM operator: vec(Zr, c->col) */
	MPI_Allreduce(MPI_IN_PLACE, Zr, c->col, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	for (i=0; i<num_of_labels; ++i) {
	  _D2_FUNC(irms)(dim, str, c->p_supp + i*strxdim, Zr + i*str);
	}

	// re-calculate C
	calculate_distmat(data_ph, label, size, c, C);
	/* rho is an important hyper-parameter */
	for (i=0; i<str*col; ++i) C[i] /= rho; // normalize C and Y

	break;

      case D2_HISTOGRAM :
	break;

      case D2_N_GRAM :
	if (iter > 0) {
	  //assert(num_of_labels * (strxdim * data_ph->vocab_size + 1) <= size);

	for (i=0; i<num_of_labels*strxdim*data_ph->vocab_size; ++i) Zr2[i] = 0; //reset Zr to temporarily storage
	for (i=0; i<size; ++i) {
	  accumulate_symbolic(dim, str, p_str[i], 
			      p_supp_sym + dim*p_str_cum[i], 
			      Z + str*p_str_cum[i], 
			      Zr2 + label[i]*strxdim*data_ph->vocab_size,
			      data_ph->vocab_size);	  
	}
#ifdef __USE_MPI__
	/* ALLREDUCE by SUM operator: vec(Zr2, num_of_labels*str*dim*data_ph->vocab_size) */
	MPI_Allreduce(MPI_IN_PLACE, Zr2, c->col * dim * data_ph->vocab_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	for (i=0; i<num_of_labels; ++i) {
	  minimize_symbolic(dim, str, c->p_supp_sym + i*strxdim, Zr2 + i*strxdim*data_ph->vocab_size, data_ph->vocab_size, data_ph->dist_mat, Zr2 + num_of_labels*strxdim*data_ph->vocab_size);
	}

	// re-calculate C
	calculate_distmat(data_ph, label, size, c, C);
	/* rho is an important hyper-parameter */
	for (i=0; i<str*col; ++i) C[i] /= rho; // normalize C and Y
	}
	break;
      }
    }    

    /*************************************************************************/
    // step 6: check residuals
    if ((iter%100==99 ) || (iter < 100 && iter%20 == 19))  {
      obj = _D2_CBLAS_FUNC(dot)(str*col, C, 1, X, 1);
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, X, 1);
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, Z0,1);
      primres = _D2_CBLAS_FUNC(asum)(str*col, X, 1);
      dualres = _D2_CBLAS_FUNC(asum)(str*col,Z0, 1);
#ifdef __USE_MPI__
      /* ALLREDUCE by SUM operator: obj, primres, dualres */
      MPI_Allreduce(MPI_IN_PLACE, &obj,     1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &primres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &dualres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
      obj     *= rho / p_data->global_size;
      primres /= p_data->global_size;
      dualres /= p_data->global_size;
      VPRINTF("\t%d\t%f\t%f\t%f\t%f\n", iter+1, obj, primres, dualres, nclock_end());
    }

    if (nclock_end() > time_budget) {break;}
  }

  _D2_FREE(Z0);
  _D2_FREE(label_count);
  if (hasZr2) _D2_FREE(Zr2);
  return 0;
}

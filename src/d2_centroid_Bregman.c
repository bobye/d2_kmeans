#include "d2_clustering.h"
#include "d2_math.h"
#include "stdio.h"

int d2_allocate_work_sphBregman(sph *ph, int size, var_sphBregman * var_phwork) {
      var_phwork->X = 
	(SCALAR *) malloc (ph->str * ph->col * sizeof(SCALAR));
      var_phwork->Z = 
	(SCALAR *) malloc (ph->str * ph->col * sizeof(SCALAR));
      var_phwork->Xc=
	(SCALAR *) malloc (ph->col * sizeof(SCALAR));
      var_phwork->Zr=
	(SCALAR *) malloc (ph->str * size * sizeof(SCALAR));
      var_phwork->Y = 
	(SCALAR *) calloc (ph->str * ph->col, sizeof(SCALAR)); // initialized
  return 0;
}

int d2_free_work_sphBregman(var_sphBregman *var_phwork) {
  free(var_phwork->X);
  free(var_phwork->Z);
  free(var_phwork->Y);
  free(var_phwork->Xc);
  free(var_phwork->Zr);
  return 0;
}

/* See matlab/centroid_sphBregman.m 
 * for a prototype implementation in Matlab. 
 */
int d2_centroid_sphBregman(mph *p_data, // data
			   var_mph * var_work,
			   int idx_ph, // index of phases
			   sph *c0,
			   __OUT__ sph *c) {
  sph *data_ph = p_data->ph + idx_ph;
  int *label = p_data->label;
  int num_of_labels = p_data->num_of_labels;
  int dim = data_ph->dim;
  int col = data_ph->col;
  int str = data_ph->str;
  int size = p_data->size;
  int *p_str = data_ph->p_str;
  SCALAR *p_supp = data_ph->p_supp;
  SCALAR *p_w = data_ph->p_w;
  SCALAR *C = var_work->g_var[idx_ph].C;
  SCALAR *X = var_work->l_var_sphBregman[idx_ph].X;
  SCALAR *Y = var_work->l_var_sphBregman[idx_ph].Y;
  SCALAR *Z = var_work->l_var_sphBregman[idx_ph].Z;
  SCALAR *Xc= var_work->l_var_sphBregman[idx_ph].Xc;
  SCALAR *Zr= var_work->l_var_sphBregman[idx_ph].Zr;

  int i,j;
  int max_niter = 1000, iter;
  SCALAR rho, obj, primres, dualres;
  SCALAR tmp, *Z0;
  SCALAR *p_scal, *p_scal2;
  int *label_count;

  /* Initialization */
  if (!c0) {
    d2_centroid_randn(p_data, idx_ph, c);// For profile purpose only
  } else {
    *c = *c0; // warm start (for clustering purpose)
  }  
  
  // compute C
  for (i=0, p_scal = C, p_scal2 = p_supp;i < size;  ++i) {
    _D2_FUNC(pdist2)(dim, str, p_str[i], &c->p_supp[label[i]*dim*str], p_scal2, p_scal);
    p_scal += str*p_str[i]; p_scal2 += dim*p_str[i];
  }

  // rho is an important hyper-parameter
  rho = _D2_CBLAS_FUNC(asum)(str*col, C, 1) / (str*col);


  /* Currently we reinitialize Z, which may not be necessary
   * Indeed, we may only need to reinitialize for entries 
   * whose label are changed  
   */
  for (i=0, p_scal = Z; i<size; ++i) {
    tmp = 1./(str * p_str[i]);
    for (j=0; j<str*p_str[i]; ++j, ++p_scal) *p_scal = tmp;
  }

  // allocate buffer of Z
  Z0 = (SCALAR *) malloc(str*col * sizeof(SCALAR));

  // calculate labels counts
  label_count = (int *) calloc(num_of_labels, sizeof(int));    
  for (i=0; i<size; ++i) 
    ++label_count[label[i]];

  // main loop
  for (iter=0; iter < max_niter; ++iter) {

    // step 1: update X    
    _D2_CBLAS_FUNC(copy)(str*col, C, 1, X, 1);
    _D2_CBLAS_FUNC(axpy)(str*col, 1, Y, 1, X, 1);
    _D2_CBLAS_FUNC(scal)(str*col, -1./rho, X, 1);
    _D2_FUNC(exp)(str*col, X);
    _D2_FUNC(vmul)(str*col, X, Z, X);
    // normalize X
    _D2_FUNC(cnorm)(str, col, X, Xc);
    _D2_FUNC(grms)(str, col, X, p_w);

    // step 2: update Z
    // Z = X.*exp(Y/rho)
    _D2_CBLAS_FUNC(copy)(str*col, Z, 1, Z0, 1);
    _D2_CBLAS_FUNC(copy)(str*col, Y, 1, Z, 1);
    _D2_CBLAS_FUNC(scal)(str*col, 1./rho, Z, 1);
    _D2_FUNC(exp)(str*col, Z);
    _D2_FUNC(vmul)(str*col, Z, X, Z);
    
    // normalize Z
    for (i=0, p_scal = Z;i<size; ++i) {
      _D2_FUNC(rnorm)(str, p_str[i], p_scal, Zr + str*i);
      _D2_FUNC(gcms)(str, p_str[i], p_scal, c->p_w + str*label[i]);
      p_scal = p_scal + str*p_str[i];
    }
   
    // step 3: update Y
    _D2_CBLAS_FUNC(axpy)(str*col, rho, X, 1, Y, 1);
    _D2_CBLAS_FUNC(axpy)(str*col, -rho, Z, 1, Y, 1);


    // step 4: update c->p_w
    for (i=0; i<str*num_of_labels; ++i) c->p_w[i] = 0; //reset c->p_w
    _D2_FUNC(cnorm)(str, size, Zr, NULL);

    for (i=0; i<size; ++i) 
      _D2_CBLAS_FUNC(axpy)(str, 1./label_count[label[i]], Zr + str*i, 1, c->p_w +label[i]*str, 1);
    _D2_FUNC(cnorm)(str, num_of_labels, c->p_w, NULL);
    

    // step 5: update c->p_supp (optional)
    for (i=0; i<str*dim*num_of_labels; ++i) c->p_supp[i] = 0; // reset c->p_supp
    for (i=0, p_scal = X, p_scal2 = p_supp;i < size;  ++i) {
      _D2_CBLAS_FUNC(gemm)(CblasColMajor, CblasNoTrans, CblasTrans, 
			   dim, str, p_str[i], 1, p_scal2, dim, p_scal, str, 1, 
			   c->p_supp + label[i]*dim*str, dim);
      p_scal += str*p_str[i]; p_scal2 += dim*p_str[i];
    }
    for (i=0; i<num_of_labels; ++i) {
      _D2_CBLAS_FUNC(scal)(str*dim, 1./label_count[i], c->p_supp + i*dim*str, 1);
      _D2_FUNC(irms)(dim, str, c->p_supp + i*dim*str, c->p_w + i*str);
    }


    // re-calculate C
    for (i=0, p_scal = C, p_scal2 = p_supp;i < size;  ++i) {
      _D2_FUNC(pdist2)(dim, str, p_str[i], c->p_supp + label[i]*dim*str, p_scal2, p_scal);
      p_scal += str*p_str[i]; p_scal2 += dim*p_str[i];
    }
    
    // step 6: check residuals
    obj = _D2_CBLAS_FUNC(dot)(str*col, C, 1, X, 1) / size;
    if (iter%100 == 0) {
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, X, 1);
      _D2_CBLAS_FUNC(axpy)(str*col, -1, Z, 1, Z0,1);
      primres = _D2_CBLAS_FUNC(asum)(str*col, X, 1) / size;
      dualres = _D2_CBLAS_FUNC(asum)(str*col,Z0, 1) / size;
      printf("\t%d\t%f\t%f\t%f\n", iter, obj, primres, dualres);
    }
  }

  free(label_count);

  return 0;
}

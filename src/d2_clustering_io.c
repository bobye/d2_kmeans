#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"
#include "d2_math.h"
#include "d2_centroid_util.h"
#ifdef __USE_MPI__
#include <mpi.h>
#endif

struct timespec io_time;

/** Load Data Set: see specification of format at README.md */
int d2_read(const char* filename, mph *p_data) {
  char filename_main[255];
  nclock_start_p(&io_time);

  if (nprocs > 1) {
    sprintf(filename_main, "%s.%d", filename, world_rank);
  } else {
    sprintf(filename_main, "%s", filename);
  }
  FILE *fp = fopen(filename_main, "r+");
  assert(fp);

  size_t i;
  int n;
  int **p_str, **p_supp_sym, str;
  double **p_supp, **p_w;
  int s_ph = p_data->s_ph;
  size_t size = p_data->size;

  p_str  = (int **) malloc(s_ph * sizeof(int *));
  p_supp_sym=(int**) malloc(s_ph * sizeof(int *));
  p_supp = (double **) malloc(s_ph * sizeof(double *));
  p_w    = (double **) malloc(s_ph * sizeof(double *));
  
  for (n=0; n<s_ph; ++n) {
    p_str[n]  = p_data->ph[n].p_str;
    p_supp_sym[n]  = p_data->ph[n].p_supp_sym;
    p_supp[n] = p_data->ph[n].p_supp;
    p_w[n]    = p_data->ph[n].p_w;
    p_data->ph[n].col = 0;
  }


  // Load header information if available
  for (n=0; n<s_ph; ++n) {
    if (p_data->ph[n].metric_type == D2_HISTOGRAM) {
      char filename_extra[255];
      FILE *fp_new; // local variable
      int str, c, i;
      sprintf(filename_extra, "%s.hist%d", filename, n);
      fp_new = fopen(filename_extra, "r+"); assert(fp_new);
      c=fscanf(fp_new, "%d", &str); assert(c>0 && str == p_data->ph[n].str);
      for (i=0; i< str*str; ++i) 
	fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data->ph[n].dist_mat[i]));
      p_data->ph[n].vocab_size = str;
      fclose(fp_new);
    }
    else if (p_data->ph[n].metric_type == D2_WORD_EMBED) {
      char filename_extra[255];
      FILE *fp_new; // local variable
      int dim, c, i;
      sprintf(filename_extra, "%s.vocab%d", filename, n);
      fp_new = fopen(filename_extra, "r+"); assert(fp_new);
      c=fscanf(fp_new, "%d", &dim); assert(c>0 && dim == p_data->ph[n].dim);
      c=fscanf(fp_new, "%d", &p_data->ph[n].vocab_size);
      p_data->ph[n].vocab_vec = _D2_MALLOC_SCALAR(dim * p_data->ph[n].vocab_size);
      for (i=0; i<p_data->ph[n].vocab_size * dim; ++i) 
	fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data->ph[n].vocab_vec[i]));
      fclose(fp_new);
    }
  }

  // Read main data file
  for (i=0; i<size; ++i) {
    for (n=0; n<s_ph; ++n) {      
      double *p_supp_sph, *p_w_sph, w_sum;
      int *p_supp_sym_sph;
      int dim, strxdim, c, j;
      // read dimension and stride    
      c=fscanf(fp, "%d", &dim); 
      if (c!=1) {
	printf("rank %d warning: only read %zd d2!\n", world_rank, i);
	p_data->size = i;
	size = i; break;
      }
      assert(dim == p_data->ph[n].dim);    
      fscanf(fp, "%d", p_str[n]); 
      str = *(p_str[n]); assert(str > 0);
      if (str >p_data->ph[n].max_str) p_data->ph[n].max_str = str;

      // check if needed to reallocate
      if (p_data->ph[n].col + str >= p_data->ph[n].max_col) {
	printf("rank %d warning: preallocated memory for phase %d is insufficient! Reallocated.\n", world_rank, n);
	if (p_data->ph[n].metric_type == D2_EUCLIDEAN_L2) {
	  p_data->ph[n].p_supp = (double *) realloc(p_data->ph[n].p_supp, 2 * dim * p_data->ph[n].max_col * sizeof(double));
	  p_data->ph[n].p_w = (double *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(double));
	  assert(p_data->ph[n].p_supp != NULL && p_data->ph[n].p_w != NULL);
	  p_supp[n] = p_data->ph[n].p_supp + p_data->ph[n].col * dim;
	  p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
	} else if (p_data->ph[n].metric_type == D2_HISTOGRAM) {
	  p_data->ph[n].p_w = (double *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(double));
	  assert(p_data->ph[n].p_w != NULL);
	  p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
	} else if (p_data->ph[n].metric_type == D2_WORD_EMBED) {
	  p_data->ph[n].p_supp_sym = (int *) realloc(p_data->ph[n].p_supp_sym, 2* p_data->ph[n].max_col * sizeof(int));
	  p_data->ph[n].p_w = (SCALAR *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(SCALAR));
	  assert(p_data->ph[n].p_supp_sym != NULL && p_data->ph[n].p_w != NULL);
	  p_supp_sym[n] = p_data->ph[n].p_supp_sym + p_data->ph[n].col;
	  p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
	}

	p_data->ph[n].max_col *= 2;		// resize
      }

      p_data->ph[n].col += str; //incr

      // read weights      
      p_w_sph = p_w[n]; w_sum = 0.;
      for (j=0; j<str; ++j) {
	fscanf(fp, SCALAR_STDIO_TYPE, &p_w_sph[j]); 
	if (p_data->ph[n].metric_type != D2_HISTOGRAM) assert(p_w_sph[j] > 1E-9);
	w_sum += p_w_sph[j];
      }
      //assert(fabs(w_sum - 1.0) <= 1E-6);
      for (j=0; j<str; ++j) {
	p_w_sph[j] /= w_sum; // re-normalize all weights
      }
      p_w[n] = p_w[n] + str;

      // read support vec
      if (p_data->ph[n].metric_type == D2_EUCLIDEAN_L2) {
	p_supp_sph = p_supp[n];strxdim = str*dim;
	for (j=0; j<strxdim; ++j)
	  fscanf(fp, SCALAR_STDIO_TYPE, &p_supp_sph[j]); 
	p_supp[n] = p_supp[n] + strxdim;
      } else if (p_data->ph[n].metric_type == D2_WORD_EMBED) {	
	p_supp_sym_sph = p_supp_sym[n]; 
	for (j=0; j<str; ++j) {
	  fscanf(fp, "%d", &p_supp_sym_sph[j]); p_supp_sym_sph[j] --; // index started at one
	  if (p_supp_sym_sph[j] < 0) p_supp_sym_sph[j] = p_data->ph[n].vocab_size - 1;
	}
	p_supp_sym[n] = p_supp_sym[n] + str;
      }
      p_str[n] ++;
    }
  }

  for (n=0; n<s_ph; ++n) {
    size_t * p_str_cum = p_data->ph[n].p_str_cum;
    int * p_str = p_data->ph[n].p_str;
    p_str_cum[0] = 0;
    for (i=1; i<size; ++i) {
      p_str_cum[i] = p_str_cum[i-1] + p_str[i-1];
    }
  }

  // free the pointer space
  free(p_w); free(p_supp); free(p_str); free(p_supp_sym);
  fclose(fp);

#ifdef __USE_MPI__
  assert(sizeof(size_t)  == sizeof(unsigned long long));
  MPI_Allreduce(&p_data->size, &p_data->global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
#else
  p_data->global_size = p_data->size;
#endif

  VPRINTF("IO time: %lf\n", nclock_end_p(&io_time));

  return 0;
}

int d2_write(const char* filename, mph *p_data) {
  FILE *fp = NULL;
  size_t i;
  const int s_ph = p_data->s_ph;
  const size_t size = p_data->size;
  if (world_rank == 0) {
  fp = filename? fopen(filename, "w+") : stdout;
  assert(fp);

  for (i=0; i<size; ++i) {
    int j;
    for (j=0; j<s_ph; ++j) 
      if (p_data->ph[j].col > 0) {
	int k, d;
	int dim = p_data->ph[j].dim;
	int str = p_data->ph[j].p_str[i];
	size_t pos = p_data->ph[j].p_str_cum[i];
	fprintf(fp, "%d\n", dim);
	fprintf(fp, "%d\n", str);
	for (k=0; k<str; ++k) fprintf(fp, "%lf ", p_data->ph[j].p_w[pos + k]);
	fprintf(fp, "\n");
	for (k=0; k<str; ++k) {
	  for (d=0; d<dim; ++d) fprintf(fp, "%lf ", p_data->ph[j].p_supp[(pos+k)*dim + d]);
	  fprintf(fp, "\n"); 
	}
      }
  }

  if (filename) fclose(fp);
  }
  return 0;
}

int d2_write_labels(const char* filename, mph *p_data) {
  FILE *fp = NULL;
  size_t i;
  int k;

  assert(filename);

  for (k=0; k<nprocs; ++k) {
    if (k == world_rank) {
      if (k==0) fp = fopen(filename, "w"); else fp=fopen(filename, "a");
      for (i=0; i<p_data->size; ++i) {
	fprintf(fp, "%d\n", p_data->label[i]);
      }
      fclose(fp);
    }
#ifdef __USE_MPI__
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  return 0;
}

/**
   Serial function to split data
 */
int d2_write_split(const char* filename, mph *p_data, int splits) {
  const int s_ph = p_data->s_ph;
  const size_t size = p_data->size; 
  size_t *indices, batch_size, n;
  int k;
  FILE *fp;
  char local_filename[255];

  
  assert(filename != NULL);

  indices = _D2_MALLOC_SIZE_T(size);
  for (n = 0; n < size; ++n) indices[n] = n; shuffle(indices, size);

  // output indices
  sprintf(local_filename, "%s.ind", filename);
  fp = fopen(local_filename, "w+"); assert(fp);
  for (n = 0; n < size; ++n) fprintf(fp, "%zd\n", indices[n]);
  fclose(fp);
  fprintf(stderr, "\twrite %zd indices to %s\n", size, local_filename);

  // output reads in several segments
  batch_size = 1 + (size-1) / splits;
  VPRINTF("batch_size: %zd\n", batch_size);

  for (k=0; k<splits; ++k) {
    size_t idx;
    FILE *fp = NULL;
    char local_filename[255];
    sprintf(local_filename, "%s.%d", filename, k);

    fp = fopen(local_filename, "w+");
    assert(fp);
    
    for (idx=batch_size*k; idx<size && idx<batch_size*(k+1); ++idx) {
      int j;
      size_t i = indices[idx];
      for (j=0; j<s_ph; ++j) 
	if (p_data->ph[j].col > 0) {
	  int k, d;
	  int dim = p_data->ph[j].dim;
	  int str = p_data->ph[j].p_str[i];
	  size_t pos = p_data->ph[j].p_str_cum[i];
	  if (1) { // non pre-processed
	    if (p_data->ph[j].metric_type == D2_HISTOGRAM) {
	      fprintf(fp, "%d %d ", dim, str);
	    } else {
	      fprintf(fp, "%d\n%d\n", dim, str);
	    }
	  for (k=0; k<str; ++k) fprintf(fp, "%lf ", p_data->ph[j].p_w[pos + k]);
	  fprintf(fp, "\n");
	  if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2) {
	    for (k=0; k<str; ++k) {
	      for (d=0; d<dim; ++d) fprintf(fp, "%lf ", p_data->ph[j].p_supp[(pos+k)*dim + d]);
	      fprintf(fp, "\n"); 
	    }
	  }
	  else if (p_data->ph[j].metric_type == D2_WORD_EMBED) {
	    for (k=0; k<str; ++k) {
	      fprintf(fp, "%d ", p_data->ph[j].p_supp_sym[pos + k]);
	    }
	    fprintf(fp, "\n");
	  }
	  } else {
	    /** Pre-processed
	     */
	    SCALAR *supp=NULL, *c_supp=NULL;
	    SCALAR *w=NULL, *c_w = NULL;
	    char renew=0;
	    if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2) {
	      supp=p_data->ph[j].p_supp + pos*dim;
	      w=p_data->ph[j].p_w + pos;
	    }
	    else if (p_data->ph[j].metric_type == D2_WORD_EMBED) {
	      renew = 1;
	      supp = (SCALAR*) _D2_MALLOC_SCALAR(str*dim);	 
	      w = (SCALAR*) _D2_MALLOC_SCALAR(str);
	      for (k=0; k<str; ++k) {
		for (d=0; d<dim; ++d) 
		  supp[k*dim + d] = p_data->ph[j].vocab_vec[p_data->ph[j].p_supp_sym[pos + k]*dim + d];
		w[k] = p_data->ph[j].p_w[pos + k];
	      }
	    }

	    if (str > p_data->ph[j].str) {	      
	      c_supp = (SCALAR*) _D2_MALLOC_SCALAR(p_data->ph[j].str*dim);
	      c_w = (SCALAR*) _D2_MALLOC_SCALAR(p_data->ph[j].str);
	      merge(dim, supp, w, str, c_supp, c_w, p_data->ph[j].str);
	      str=p_data->ph[j].str;
	    } else {
	      c_supp = supp;
	      c_w = w;
	    }

	    fprintf(fp, "%d\n", dim);
	    fprintf(fp, "%d\n", str);
	    for (k=0; k<str; ++k) fprintf(fp, "%lf ", c_w[k]);
	    fprintf(fp, "\n");	  
	    for (k=0; k<str; ++k) {
	      for (d=0; d<dim; ++d) fprintf(fp, "%lf ", c_supp[k*dim + d]);
	      fprintf(fp, "\n"); 
	    }

	    
	    if (renew) _D2_FREE(supp); 
	    if (renew) _D2_FREE(w); 
	    if (c_supp!=supp) _D2_FREE(c_supp);
	    if (c_w!=w) _D2_FREE(c_w);
	  }
	}
    }

    fclose(fp);
    fprintf(stderr, "\twrite %zd objects to %s\n", idx - batch_size * k, local_filename);
  }
  return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "d2/clustering.h"
#include "d2/math.h"
#include "d2/centroid_util.h"
#ifdef __USE_MPI__
#include <mpi.h>
#endif

/** Load Data Set: see specification of format at README.md */
int d2_read(const char* filename, const char* meta_filename, mph *p_data) {
  char filename_main[255];
  FILE *fp =NULL;
  double io_startTime;
  io_startTime = getRealTime();

  if (nprocs > 1) {
    sprintf(filename_main, "%s.%d", filename, world_rank);
    fp = fopen(filename_main, "r");
  }

  if (nprocs == 1 || !fp ) {
    sprintf(filename_main, "%s", filename);
    fp = fopen(filename_main, "r");
  }

  assert(fp);

  size_t i;
  int n;
  int **p_str, **p_supp_sym;
  SCALAR **p_supp, **p_w;
  int s_ph = p_data->s_ph;
  size_t size = p_data->size;

  p_str  = (int **) malloc(s_ph * sizeof(int *));
  p_supp_sym=(int**) malloc(s_ph * sizeof(int *));
  p_supp = (SCALAR **) malloc(s_ph * sizeof(SCALAR *));
  p_w    = (SCALAR **) malloc(s_ph * sizeof(SCALAR *));
  
  for (n=0; n<s_ph; ++n) {
    p_str[n]  = p_data->ph[n].p_str;
    p_supp_sym[n]  = p_data->ph[n].p_supp_sym;
    p_supp[n] = p_data->ph[n].p_supp;
    p_w[n]    = p_data->ph[n].p_w;
    p_data->ph[n].col = 0;
  }


  // Load header information if available
  for (n=0; n<s_ph; ++n) {
    if (p_data->ph[n].metric_type == D2_HISTOGRAM ||
	p_data->ph[n].metric_type == D2_SPARSE_HISTOGRAM) {
      char filename_extra[255];
      FILE *fp_new; // local variable
      int str, c, i;
      if (meta_filename && s_ph == 1) 
	strcpy(filename_extra, meta_filename);
      else
	sprintf(filename_extra, "%s.hist%d", filename, n);
      fp_new = fopen(filename_extra, "r"); assert(fp_new);
      c=fscanf(fp_new, "%d", &str); //assert(c>0 && str == p_data->ph[n].str);
      p_data->ph[n].dist_mat = _D2_MALLOC_SCALAR(str * str);
      p_data->ph[n].is_meta_allocated = true;
      for (i=0; i< str*str; ++i) 
	fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data->ph[n].dist_mat[i]));
      p_data->ph[n].vocab_size = str;
      fclose(fp_new);
    }
    else if (p_data->ph[n].metric_type == D2_WORD_EMBED) {
      char filename_extra[255];
      FILE *fp_new; // local variable
      int dim, c, i;
      if (meta_filename && s_ph == 1) 
	strcpy(filename_extra, meta_filename);
      else
	sprintf(filename_extra, "%s.vocab%d", filename, n);
      fp_new = fopen(filename_extra, "r"); assert(fp_new);
      c=fscanf(fp_new, "%d", &dim); assert(c>0 && dim == p_data->ph[n].dim);
      c=fscanf(fp_new, "%d", &p_data->ph[n].vocab_size);
      p_data->ph[n].vocab_vec = _D2_MALLOC_SCALAR(dim * p_data->ph[n].vocab_size);
      p_data->ph[n].is_meta_allocated = true;
      for (i=0; i<p_data->ph[n].vocab_size * dim; ++i) 
	fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data->ph[n].vocab_vec[i]));
      fclose(fp_new);
    }
  }

  // Read main data file
  for (i=0; i<size; ++i) {
    for (n=0; n<s_ph; ++n) {      
      SCALAR *p_supp_sph, *p_w_sph, w_sum;
      int *p_supp_sym_sph;
      int dim, str, strxdim, c, j;
      // read dimension and stride    
      c=fscanf(fp, "%d", &dim); 
      if (c!=1) {
	printf("rank %d warning: only read %zd d2!\n", world_rank, i);
	p_data->size = i;
	size = i; break;
      }
      assert(dim == p_data->ph[n].dim);    
      fscanf(fp, "%d", &str); assert(str >= 0);
      if (str == 0) continue;
      *p_str[n] = str; 

      if (str >p_data->ph[n].max_str) p_data->ph[n].max_str = str;

      // check if needed to reallocate
      if (p_data->ph[n].col + str > p_data->ph[n].max_col) {
	printf("rank %d warning: preallocated memory for phase %d is insufficient! Reallocated.\n", world_rank, n);
	if (p_data->ph[n].metric_type == D2_EUCLIDEAN_L2) {
	  p_data->ph[n].p_supp = (SCALAR *) realloc(p_data->ph[n].p_supp, 2 * dim * p_data->ph[n].max_col * sizeof(SCALAR));
	  p_data->ph[n].p_w = (SCALAR *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(SCALAR));
	  assert(p_data->ph[n].p_supp != NULL && p_data->ph[n].p_w != NULL);
	  p_supp[n] = p_data->ph[n].p_supp + p_data->ph[n].col * dim;
	  p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
	} else if (p_data->ph[n].metric_type == D2_HISTOGRAM) {
	  p_data->ph[n].p_w = (SCALAR *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(SCALAR));
	  assert(p_data->ph[n].p_w != NULL);
	  p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
	} else if (p_data->ph[n].metric_type == D2_WORD_EMBED ||
		   p_data->ph[n].metric_type == D2_SPARSE_HISTOGRAM ) {
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
	/* Unless the input is of format D2_HISTOGRAM,
	   we require all support points should have non-zero weights */
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
      } else if (p_data->ph[n].metric_type == D2_WORD_EMBED ||
		 p_data->ph[n].metric_type == D2_SPARSE_HISTOGRAM) {	
	p_supp_sym_sph = p_supp_sym[n]; 
	for (j=0; j<str; ++j) {
	  fscanf(fp, "%d", &p_supp_sym_sph[j]); p_supp_sym_sph[j] --; // index read started at one
	  if (p_supp_sym_sph[j] < 0)
	    p_supp_sym_sph[j] = p_data->ph[n].vocab_size - 1; // handle boundary case
	}
	p_supp_sym[n] = p_supp_sym[n] + str;
      }
      p_str[n] ++;
    }
  }

  for (n=0; n<s_ph; ++n) 
  if (p_data->ph[n].col > 0) {
    size_t * p_str_cum = p_data->ph[n].p_str_cum;
    int * p_str = p_data->ph[n].p_str;
    p_str_cum[0] = 0;
    for (i=1; i<size; ++i) {
      p_str_cum[i] = p_str_cum[i-1] + p_str[i-1];
    }
    if (p_data->ph[n].metric_type == D2_SPARSE_HISTOGRAM) {
      // change str to vocab_size, which is for initializing centroids.
      p_data->ph[n].str = p_data->ph[n].vocab_size; 
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

  VPRINTF("IO time: %lf\n", getRealTime() - io_startTime);

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
	if (p_data->ph[j].metric_type == D2_HISTOGRAM) continue;
	if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2) {
	  for (k=0; k<str; ++k) {
	    for (d=0; d<dim; ++d) fprintf(fp, "%lf ", p_data->ph[j].p_supp[(pos+k)*dim + d]);
	    fprintf(fp, "\n"); 
	  }
	}
	if (p_data->ph[j].metric_type == D2_WORD_EMBED ||
	    p_data->ph[j].metric_type == D2_SPARSE_HISTOGRAM) {
	  for (k=0; k<str; ++k) {
	    fprintf(fp, "%d ", p_data->ph[j].p_supp_sym[pos + k] + 1);	    
	  }
	  fprintf(fp, "\n");
	}
      } else {
	fprintf(fp, "%d\n", p_data->ph[j].dim);
	fprintf(fp, "0\n");
	fprintf(fp, "\n");
      }
  }

  if (filename) fclose(fp);
  VPRINTF("Write centroids to %s\n", filename);
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

  VPRINTF("Write data partitioned labels to %s\n", filename);

  return 0;
}

int d2_write_labels_serial(const char* filename_ind, const char* filename, mph *p_data) {
  FILE *fp = NULL;
  size_t i, global_size = p_data->global_size;

  assert(filename);
  
  if (0 == world_rank) {
    int *label, *label_o;
    size_t *indice;
    char filename_label[255];
    label = _D2_MALLOC_INT(2*global_size); label_o = label + global_size;
    indice = _D2_MALLOC_SIZE_T(global_size);

    sprintf(filename_label, "%s.label", filename);
    fp = fopen(filename_label, "r"); assert(fp);
    for (i=0; i<global_size; ++i) {
      fscanf(fp, "%d\n", &label[i]);
    }
    fclose(fp);

    fp = fopen(filename_ind, "r"); 
    if (fp) { // if ind file exists => has data partition 
      for (i=0; i<global_size; ++i) {
	fscanf(fp, "%zu\n", &indice[i]);
      }
      fclose(fp);

      for (i=0; i<global_size; ++i) {
	label_o[indice[i]] = label[i];
      }

      sprintf(filename_label, "%s.label_o", filename);
      fp = fopen(filename_label, "w");
      for (i=0; i<global_size; ++i) {
	fprintf(fp, "%d\n", label_o[i]);
      }
      fclose(fp);
      VPRINTF("Write serial data labels to %s\n", filename_label);
    }
    _D2_FREE(label); _D2_FREE(indice);
  }


#ifdef __USE_MPI__
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  return 0;  
}

/**
   Serial function to split data
 */
int d2_write_split(const char* filename, mph *p_data, int splits, char is_pre_processed) {
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
	  if (!is_pre_processed) { // non pre-processed
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
	    else if (p_data->ph[j].metric_type == D2_WORD_EMBED ||
		     p_data->ph[j].metric_type == D2_SPARSE_HISTOGRAM) {
	      for (k=0; k<str; ++k) {
		fprintf(fp, "%d ", p_data->ph[j].p_supp_sym[pos + k] + 1); // fix an important index bug
	      }
	      fprintf(fp, "\n");
	    }
	  } else {
	    /** Pre-processed
	     */
	    SCALAR *supp=NULL, *w=NULL;
	    int *supp_sym=NULL;
	    SCALAR *p_supp=NULL, *p_w=NULL;
	    int *p_supp_sym=NULL;


	    w=p_data->ph[j].p_w + pos;
	    if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2) {
	      supp=p_data->ph[j].p_supp + pos*dim;
	      p_supp = supp;
	      p_w = w;
	    }	    

	    if (p_data->ph[j].metric_type == D2_WORD_EMBED) {
	      supp_sym = p_data->ph[j].p_supp_sym + pos;
	      p_supp = (SCALAR*) _D2_MALLOC_SCALAR(str*dim);	 
	      p_w = (SCALAR*) _D2_MALLOC_SCALAR(str);
	      for (k=0; k<str; ++k) {
		for (d=0; d<dim; ++d) 
		  p_supp[k*dim + d] = p_data->ph[j].vocab_vec[supp_sym[k]*dim + d];
		p_w[k] = w[k];
	      }

	    }

	    if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2 ||
		p_data->ph[j].metric_type == D2_WORD_EMBED) {
	    if (str > p_data->ph[j].str) {
	      // This preprocessing makes a D2 into one with support points less than p_data->ph[j].str;
	      SCALAR *c_supp = (SCALAR*) _D2_MALLOC_SCALAR(p_data->ph[j].str*dim);
	      SCALAR *c_w = (SCALAR*) _D2_MALLOC_SCALAR(p_data->ph[j].str);
	      merge(dim, p_supp, p_w, str, c_supp, c_w, p_data->ph[j].str);
	      str=p_data->ph[j].str;
	      if (p_supp!=supp) _D2_FREE(p_supp);
	      if (p_w!=w) _D2_FREE(p_w);
	      p_supp = c_supp;
	      p_w = c_w;
	    }
	    }

	    if (p_data->ph[j].metric_type == D2_HISTOGRAM) {
	      // This preprocessing makes a dense histogram into a sparse one
	      int nnz = 0, ind = 0;
	      for (k=0; k<str; ++k) nnz += (w[k] > 0);
	      p_supp_sym = (int*) _D2_MALLOC_INT(nnz);
	      p_w = (SCALAR*) _D2_MALLOC_SCALAR(nnz);
	      for (k=0; k<str; ++k)
		if (w[k] > 0) {		  
		  p_supp_sym[ind] = k;
		  p_w[ind] = w[k];
		  ++ind;
		}
	      str = nnz; // change str to sparse counts
	    }


	    fprintf(fp, "%d\n", dim);
	    fprintf(fp, "%d\n", str);
	    for (k=0; k<str; ++k) fprintf(fp, "%lf ", p_w[k]);
	    fprintf(fp, "\n");
	    if (p_data->ph[j].metric_type == D2_EUCLIDEAN_L2 ||
		p_data->ph[j].metric_type == D2_WORD_EMBED) {
	      for (k=0; k<str; ++k) {
		for (d=0; d<dim; ++d) fprintf(fp, "%lf ", p_supp[k*dim + d]);
		fprintf(fp, "\n"); 
	      }
	    } else if (p_data->ph[j].metric_type == D2_HISTOGRAM) {
	      for (k=0; k<str; ++k) {
		fprintf(fp, "%d ", p_supp_sym[k] + 1);
	      }
	      fprintf(fp, "\n");
	    }
	    
	    if (p_supp!=supp && p_supp) _D2_FREE(p_supp);
	    if (p_w!=w && p_w) _D2_FREE(p_w);
	    if (p_supp_sym != supp_sym && p_supp_sym) _D2_FREE(p_supp_sym);
	  }
	}
    }

    fclose(fp);
    fprintf(stderr, "\twrite %zd objects to %s\n", idx - batch_size * k, local_filename);
  }
  return 0;
}





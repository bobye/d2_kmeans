#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "d2_clustering.h"
#include "d2_math.h"
#ifdef __USE_MPI__
#include <mpi.h>
#endif

struct timespec io_time;

/** Load Data Set: see specification of format at README.md */
int d2_read(const char* filename, mph *p_data) {
  nclock_start_p(&io_time);

  FILE *fp = fopen(filename, "r+");
  assert(fp);

  size_t i;
  int n;
  int **p_str, str;
  double **p_supp, **p_w;
  int s_ph = p_data->s_ph;
  size_t size = p_data->size;

  p_str  = (int **) malloc(s_ph * sizeof(int *));
  p_supp = (double **) malloc(s_ph * sizeof(double *));
  p_w    = (double **) malloc(s_ph * sizeof(double *));
  
  for (n=0; n<s_ph; ++n) {
    p_str[n]  = p_data->ph[n].p_str;
    p_supp[n] = p_data->ph[n].p_supp;
    p_w[n]    = p_data->ph[n].p_w;
    p_data->ph[n].col = 0;
  }


  for (n=0; n<s_ph; ++n) {
    if (p_data->ph[n].dim == 0) {
      char filename_extra[255];
      FILE *fp_new; // local variable
      int str, c, i;
      sprintf(filename_extra, "%s.hist%d", filename, n);
      fp_new = fopen(filename_extra, "r+"); assert(fp_new);
      c=fscanf(fp_new, "%d", &str); assert(c>0 && str == p_data->ph[n].str);
      for (i=0; i< str*str; ++i) fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data->ph[n].dist_mat[i]));
      fclose(fp_new);
    }
  }

  
  for (i=0; i<size; ++i) {
    for (n=0; n<s_ph; ++n) {      
      double *p_supp_sph, *p_w_sph, w_sum;
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

      if (p_data->ph[n].col + str >= p_data->ph[n].max_col) {
	VPRINTF("Warning: preallocated memory for phase %d is insufficient! Reallocated.\n", n);
	p_data->ph[n].p_supp = (double *) realloc(p_data->ph[n].p_supp, 2 * dim * p_data->ph[n].max_col * sizeof(double));
	p_data->ph[n].p_w = (double *) realloc(p_data->ph[n].p_w, 2* p_data->ph[n].max_col * sizeof(double));
	assert(p_data->ph[n].p_supp != NULL && p_data->ph[n].p_w != NULL);
	p_data->ph[n].max_col *= 2;		// resize
	if (dim > 0) p_supp[n] = p_data->ph[n].p_supp + p_data->ph[n].col * dim;
	p_w[n]    = p_data->ph[n].p_w + p_data->ph[n].col;
      }

      p_data->ph[n].col += str; 

      // read weights      
      p_w_sph = p_w[n]; w_sum = 0.;
      for (j=0; j<str; ++j) {
	fscanf(fp, SCALAR_STDIO_TYPE, &p_w_sph[j]); assert(p_w_sph[j] > 1E-9);
	w_sum += p_w_sph[j];
      }
      //assert(fabs(w_sum - 1.0) <= 1E-6);
      for (j=0; j<str; ++j) {
	p_w_sph[j] /= w_sum; // re-normalize all weights
      }
      p_w[n] = p_w[n] + str;

      // read support vec
      if (dim > 0) {
	p_supp_sph = p_supp[n];strxdim = str*dim;
	for (j=0; j<strxdim; ++j)
	  fscanf(fp, SCALAR_STDIO_TYPE, &p_supp_sph[j]); 
	p_supp[n] = p_supp[n] + strxdim;
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
  free(p_w); free(p_supp); free(p_str); 
  fclose(fp);

#ifdef __USE_MPI__
  assert(sizeof(size_t)  == sizeof(uint64_t));
  MPI_Allreduce(&p_data->size, &p_data->global_size, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
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


int d2_write_split(const char* filename, mph *p_data, int splits) {
  const int s_ph = p_data->s_ph;
  const size_t size = p_data->size; 
  size_t *indices, batch_size, n;
  int k;
  
  assert(filename != NULL);

  indices = _D2_MALLOC_SIZE_T(size);
  for (n = 0; n < size; ++n) indices[n] = n; shuffle(indices, size);
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

    fclose(fp);
    VPRINTF("\twrite %zd objects to %s\n", idx - batch_size * k, local_filename);
  }
  return 0;
}

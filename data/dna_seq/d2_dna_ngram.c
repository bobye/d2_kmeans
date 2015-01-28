/**
 * This file includes IO codes for reading DNA ngram data
 * and merge them with existing clustering framework. 
 *
 * @param(p_data->ph[n].p_supp_sym) The symbolic arrays.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <float.h>
#include "d2_clustering.h"
#include "d2_math.h"


/* centroid methods
 * 0: Bregman ADMM
 * 1: ADMM
 * 2: Gradient Decent
 */
#include "d2_param.h"
int d2_alg_type = D2_CENTROID_BADMM;

int key[255]; // map ASCII to index

int d2_allocate_sph_dna(sph *p_data_sph,
		    const int d,
		    const int stride,
		    const size_t num,
		    double semicol) {

  int n, m;
  assert(stride >0 && num >0);

  n = num * (stride + semicol) * d; // pre-allocate
  m = num * (stride + semicol); p_data_sph->max_col = m;

  p_data_sph->dim = d;  
  p_data_sph->str = stride;
  p_data_sph->col = 0;
  //  p_data_sph->size = num;


  p_data_sph->p_str  = _D2_CALLOC_INT(num);
  p_data_sph->p_str_cum  = _D2_CALLOC_SIZE_T(num);
  p_data_sph->p_w    = _D2_MALLOC_SCALAR(m);


  p_data_sph->p_supp_sym = _D2_MALLOC_INT(n);

  p_data_sph->vocab_size = 20; // 20 types of amino acids
  p_data_sph->dist_mat = _D2_MALLOC_SCALAR( p_data_sph->vocab_size * p_data_sph->vocab_size);


  return 0;
}



int d2_read_sph_dna(const char* filename, sph *p_data_sph, mph *p_data) {
  FILE *fp;

  size_t i,j, count_w, count_supp;
  int n;
  int s_ph = p_data->s_ph;
  size_t size = p_data->size;

  FILE *fp_new; // local variable
  int str, c;
  char symbol, name[1000];

  
  /* Begin: load keymap and their distances */
  fp_new = fopen("PAM250_distmat.dat", "r+"); assert(fp_new);
  
  n = 0; while ((c=fscanf(fp_new, "%c", &symbol)) > 0 && symbol !='\n') { 
    if (symbol >= 'A' && symbol <= 'Z') {key[symbol] = n; ++n;}
  }
  assert(n==p_data_sph->vocab_size);

  for (i=0; i< n*n; ++i) 
    { c=fscanf(fp_new, SCALAR_STDIO_TYPE, &(p_data_sph->dist_mat[i])); assert(c>0);}
  fclose(fp_new);
  /* End: PAM250_distmat.dat */

  /* Being: Load n-gram data */
  fp = fopen(filename, "r+");assert(fp);

  c=fscanf(fp, "Seq #%zd\n", &i); assert(c>=1 && i >= size);
  c=fscanf(fp, "Class #%d\n", &n); assert(c>=1);

  count_w = 0; count_supp = 0;
  for (i=0; i<size; ++i) {
    int dim, str;
    SCALAR tmp_val;
    c = fscanf(fp, "%s %d", name, &n); assert(c==2);
    c = fscanf(fp, "%d %d", &dim, &str); assert(c==2 && dim == p_data_sph->dim && str > 0);

    if (p_data_sph->col + str >= p_data_sph->max_col) {
      VPRINTF(("Warning: preallocated memory as it is insufficient! Reallocated.\n"));
      p_data_sph->p_supp_sym = (int *) realloc(p_data_sph->p_supp_sym, 2 * dim * p_data_sph->max_col * sizeof(int));
      p_data_sph->p_w = (SCALAR *) realloc(p_data_sph->p_w, 2* p_data_sph->max_col * sizeof(SCALAR));
      assert(p_data_sph->p_supp_sym != NULL && p_data_sph->p_w != NULL);
      p_data_sph->max_col *= 2;		// resize
    }

    p_data_sph->p_str_cum[i] = count_w;


    for (n=0; n<str; ++n) {
      c=fscanf(fp, SCALAR_STDIO_TYPE, &tmp_val); 
      p_data_sph->p_w[count_w++] = tmp_val;
    }


    for (n=0; n<str; ++n) {
      int m;
      for (m=0; m<dim; ++m) {
	symbol=fgetc(fp);
	if (symbol == '-') {
	  symbol = fgetc(fp); count_w --; 
	  break;
	}
	p_data_sph->p_supp_sym[count_supp*dim + m] = key[symbol];	
      }
      symbol = fgetc(fp);
      count_supp++;
    }
    
    str = str - count_supp + count_w;
    count_supp = count_w;

    p_data_sph->p_str[i] = str;
    p_data_sph->col += str;

  }

  // free the pointer space
  fclose(fp);
  return 0;
}





int d2_read_dna(mph *p_data,
		const int size_of_phases,
		const size_t size_of_samples,
		const int *avg_strides, /**
					   It is very important to make sure 
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
  //for (i=0; i<p_data->size; ++i)  p_data->label[i] = -1;

  for (i=0; i<p_data->s_ph; ++i) {
    char filename[255];
    sprintf(filename, "protein_%dgram.dat", i+1);
    printf("Load %s ...\n", filename);
    d2_allocate_sph_dna(p_data->ph + i, 
			dimension_of_phases[i], 
			avg_strides[i], 
			size_of_samples, 
			0.6);
    d2_read_sph_dna(filename, p_data->ph + i, p_data);
  }

  return success;
}



int main(int argc, char *argv[]) {
  int avg_strides[3] = {20, 20, 20};
  int dimension_of_phases[3] = {1, 2, 3};
  int size_of_phases = 3;
  size_t size_of_samples = 10742;
  int number_of_clusters = 1;
  int selected_phase = -1;
  char use_triangle = false;

  mph data;
  mph c;

  int err = d2_read_dna(&data, 
			size_of_phases,
			size_of_samples,
			&avg_strides[0],
			&dimension_of_phases[0]);

  

  d2_clustering(number_of_clusters, 
		100, 
		&data, 
		&c, 
		selected_phase,
		use_triangle);
    
  return 0;
}

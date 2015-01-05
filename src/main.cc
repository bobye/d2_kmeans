#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>

#include "d2_clustering.h"
#include "util.hh"
#include <cstdlib>
#include <getopt.h> /* for getopt_long; GNU extension */

/* centroid methods
 * 0: Bregman ADMM
 * 1: ADMM
 * 2: Gradient Decent
 */
#include "d2_param.h"
int d2_alg_type = D2_CENTROID_GRADDEC;
//int d2_alg_type = D2_CENTROID_BADMM;

int main(int argc, char *argv[])
{ 
  using namespace std;

  int size_of_phases;
  long size_of_samples;
  char *ss1_c_str = 0, *ss2_c_str = 0, *filename = 0, *ofilename = 0;

  /* default settings */
  int selected_phase = -1; 
  int number_of_clusters = 3; 
  int max_iters = 20; 

  /* IO specification */
  int ch;
  static struct option long_options[] = {
    {"strides", 1, 0, 's'},
    {"phase", 1, 0, 'p'},
    {"ifile", 1, 0, 'i'},
    {"ofile", 1, 0, 'o'},
    {"phase_only", 1, 0, 't'},
    {"clusters", 1, 0, 'c'},
    {"max_iters", 1, 0, 'm'},
    {NULL, 0, NULL, 0}
  };
  
  int option_index = 0;
  while ( (ch = getopt_long(argc, argv, "p:n:d:s:i:t:c:m:", long_options, &option_index)) != -1) {
    switch (ch) {
    case 'i':
      filename = optarg;
      break;
    case 'o':
      ofilename = optarg;
      break;
    case 'p':
      size_of_phases = atoi(optarg);
      break;
    case 'n':
      size_of_samples = atol(optarg);
      break;
    case 'd':
      ss1_c_str = optarg;
      break;
    case 's':
      ss2_c_str = optarg;
      break;
    case 't':
      selected_phase = atoi(optarg); assert(selected_phase >= 0);
      break;
    case 'c':
      number_of_clusters = atoi(optarg); assert(number_of_clusters > 0);
      break;
    case 'm':
      max_iters = max(atoi(optarg), max_iters);
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", ch);
    }
  }
  

  vector<int> dimension_of_phases(size_of_phases);  
  vector<int> avg_strides(size_of_phases);

  vector<string> ss1 = split(string(ss1_c_str), ',');
  vector<string> ss2 = split(string(ss2_c_str), ',');

  assert(size_of_phases == (int) ss1.size() && size_of_phases == (int) ss2.size());

  for (int i=0; i<size_of_phases; ++i) {
    dimension_of_phases[i] = atoi(ss1[i].c_str());
    avg_strides[i] = atoi(ss2[i].c_str());
  }    
  

  mph data;

  FILE *fp;
  
  int err = d2_allocate(&data, 
			size_of_phases,
			size_of_samples,
			&avg_strides[0],
			&dimension_of_phases[0]);

  if (err == 0) {
    fp = fopen(filename, "r+");
    d2_read(fp, &data);  
    fclose(fp);
  } else {
    cerr << "Allocation Failed!" << endl;
  }

  mph c;

  d2_clustering(number_of_clusters, max_iters, &data, &c, selected_phase);

  if (ofilename) {
    fp = fopen(ofilename, "w");
    d2_write(fp, &c);
    fclose(fp);
  } else {
    d2_write(stdout, &c);
  }
  
  d2_free(&data);
  d2_free(&c);

  return 0;
}

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
int d2_alg_type = D2_CENTROID_BADMM;

int main(int argc, char *argv[])
{ 
  using namespace std;

  int size_of_phases = 1;
  long size_of_samples;
  char *ss1_c_str = 0, *ss2_c_str = 0, *filename = 0, *ofilename = 0;
  char use_triangle = true;
  /* default settings */
  int selected_phase = -1; 
  int number_of_clusters = 3; 
  int max_iters = 50; 

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
    {"centroid_method", 1, 0, 'M'},
    {"non_triangle", 0, 0, 'T'},
    {NULL, 0, NULL, 0}
  };

  /* [BEGIN] Parsing program arguments */
  int option_index = 0;
  while ( (ch = getopt_long(argc, argv, "p:n:d:s:i:o:t:c:m:M:T", long_options, &option_index)) != -1) {
    switch (ch) {
    case 'i': /* input filename */
      filename = optarg;
      break;
    case 'o': /* output filename */
      ofilename = optarg;
      break;
    case 'p': 
      size_of_phases = atoi(optarg);
      break;
    case 'n': /* size of samples expected to be loaded */
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
    case 'M':
      d2_alg_type = atoi(optarg);
      assert(d2_alg_type == D2_CENTROID_BADMM || d2_alg_type == D2_CENTROID_GRADDEC || d2_alg_type == D2_CENTROID_ADMM);
      break;
    case 'T':
      use_triangle = false;
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", ch);
    }
  }
  

  vector<int> dimension_of_phases(size_of_phases, 0);  
  vector<int> avg_strides(size_of_phases, 0);

  vector<string> ss1 = ss1_c_str? split(string(ss1_c_str), ',') : vector<string> (1, "0");
  vector<string> ss2 = split(string(ss2_c_str), ',');

  assert(size_of_phases == (int) ss1.size() 
	 && size_of_phases == (int) ss2.size()
	 && ss2_c_str);

  cout << "Task: " << endl;
  for (int i=0; i<size_of_phases; ++i) {
    dimension_of_phases[i] = atoi(ss1[i].c_str());
    avg_strides[i] = atoi(ss2[i].c_str());
    if (dimension_of_phases[i] == 0) {
      cout << "\t" << i << "-th phase is of histogram format" << endl;
    } else if (dimension_of_phases[i] > 0) {
      cout << "\t" << i << "-th phase is of discrete distribution format" << endl;
    }
    assert(dimension_of_phases[i] >= 0 && avg_strides[i] > 0);
  }     
  /* [END] Parsing program arguments */


  /**********************************************************************************/
  /* [BEGIN] Start main program */

  /* data structure storing all information about multi-phase discrete distributions */
  mph data;
  
  int err = d2_allocate(&data, 
			size_of_phases,
			size_of_samples,
			&avg_strides[0],
			&dimension_of_phases[0]);


  if (err == 0) {
    d2_read(filename, &data);  
  } else {
    cerr << "Allocation Failed!" << endl;
  }
  
  /* data structure storing information about centroids of clusters */
  mph c; 
  c.ph = NULL; // make sure c is (re-)initialized 

  if (selected_phase >= 0 && size_of_phases > 1) {
    cout << "Clustering upon " << selected_phase <<"-th phase" << endl;
  } else if (selected_phase < 0 && size_of_phases == 1) {
    cout << "Clustering upon all phases (more than one)" << endl;
  }

  d2_clustering(number_of_clusters, 
		max_iters, 
		&data, 
		&c, 
		selected_phase,
		use_triangle);

  if (ofilename) {
    cout << "Write computed centroids to " << ofilename << endl;
    d2_write(ofilename, &c);
  } else {
    d2_write(NULL, &c);
  }
  
  d2_free(&data);
  d2_free(&c);

  cout << "[Finish!]" <<endl;
  return 0;
}

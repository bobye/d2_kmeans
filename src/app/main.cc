#include <stdio.h>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <assert.h>
#include <time.h>       /* time_t, struct tm, time, localtime */

#include "d2/clustering.h"
#include "util.hh"
#include <cstdlib>
#include <getopt.h> /* for getopt_long; GNU extension */

#ifdef __USE_MPI__
#include <mpi.h>
#endif

/* centroid methods
 * 0: Bregman ADMM
 * 1: ADMM
 * 2: Gradient Decent
 */
#include "d2/param.h"
extern int d2_alg_type;

int main(int argc, char *argv[])
{ 
#ifdef __USE_MPI__
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

  using namespace std;

  int size_of_phases = 1;
  long size_of_samples;
  char *ss1_c_str = 0, *ss2_c_str = 0, *ss3_c_str = 0,
    *filename = 0, *centroid_filename = 0, 
    *meta_filename = 0,
    *output_filename = 0,
    is_eval=0, is_load=0,
    is_pre_processed=0;
  char use_triangle = true;
  /* default settings */
  int selected_phase = -1; 
  int number_of_clusters = 3; 
  int max_iters = 100; 
  size_t num_of_batches = 0; // default not used, for prepare data only

  /* IO specification */
  int ch;
  static struct option long_options[] = {
    {"strides", 1, 0, 's'},
    {"phase", 1, 0, 'p'},
    {"ifile", 1, 0, 'i'},
    {"ofile", 1, 0, 'o'},
    {"metafile", 1, 0, 'D'},
    {"phase_only", 1, 0, 't'},
    {"clusters", 1, 0, 'k'},
    {"max_iters", 1, 0, 'm'},
    {"centroid_method", 1, 0, 'M'},
    {"non_triangle", 0, 0, 'T'},
    {"prepare_batches", 1, 0, 'P'},
    {"pre_process", 0, 0, 'Q'},
    {"types", 1, 0, 'E'},
    {"eval", 1, 0, 'e'},
    {"load", 1, 0, 'L'},
    {NULL, 0, NULL, 0}
  };

  /* [BEGIN] Parsing program arguments */
  int option_index = 0;
  while ( (ch = getopt_long(argc, argv, "p:n:s:i:o:D:d:t:k:m:M:TQP:E:e:L:", long_options, &option_index)) != -1) {
    switch (ch) {
    case 'i': /* input filename */
      filename = optarg;
      break;
    case 'D':
      meta_filename = optarg;
      break;
    case 'o':
      output_filename = optarg;
      break;
    case 'e': 
      centroid_filename = optarg;
      is_eval = 1; assert(!is_load);      
      break;
    case 'L':
      centroid_filename = optarg;
      is_load = 1; assert(!is_eval);
      break;
    case 'p': 
      size_of_phases = atoi(optarg); assert(size_of_phases > 0);
      break;
    case 'n': /* size of samples expected to be loaded */
      size_of_samples = atol(optarg); assert(size_of_samples > 0);
      break;
    case 'd': 
      ss1_c_str = optarg;
      break;
    case 's': 
      ss2_c_str = optarg;
      break;
    case 'E':
      ss3_c_str = optarg;
      break;
    case 't':
      selected_phase = atoi(optarg); assert(selected_phase >= 0);
      break;
    case 'k': 
      number_of_clusters = atoi(optarg); assert(number_of_clusters > 0);
      break;
    case 'm':
      max_iters = atoi(optarg); //max(atoi(optarg), max_iters);
      break;
    case 'M':
      d2_alg_type = atoi(optarg);
      assert(d2_alg_type == D2_CENTROID_BADMM || d2_alg_type == D2_CENTROID_GRADDEC || d2_alg_type == D2_CENTROID_ADMM);
      break;
    case 'T':
      use_triangle = false;
      break;
    case 'Q':
      is_pre_processed = true;
      break;
    case 'P':
      num_of_batches = atoi(optarg); assert(num_of_batches > 0);
      break;
    default:
      printf ("?? getopt returned character code 0%o ??\n", ch);
      exit(0);
    }
  }
  

  vector<int> dimension_of_phases(size_of_phases, 0);  
  vector<int> avg_strides(size_of_phases, 0);
  vector<int> type_of_phases(size_of_phases, 0);

  vector<string> ss1 = ss1_c_str? split(string(ss1_c_str), ',') : vector<string> (size_of_phases, "0");
  vector<string> ss2 = split(string(ss2_c_str), ',');
  vector<string> ss3 = ss3_c_str? split(string(ss3_c_str), ',') : vector<string> (size_of_phases, "0"); // default is D2_EUCLIDEAN_L2

  assert(size_of_phases == (int) ss1.size() 
	 && size_of_phases == (int) ss2.size()
	 && size_of_phases == (int) ss3.size()
	 && ss2_c_str);
  assert((size_of_phases == 1 || !meta_filename));

  if (world_rank == 0) {cout << "Task: " << endl;}  
  for (int i=0; i<size_of_phases; ++i) {
    dimension_of_phases[i] = atoi(ss1[i].c_str());
    avg_strides[i] = atoi(ss2[i].c_str());
    type_of_phases[i] = atoi(ss3[i].c_str());
    if (world_rank == 0) {
      if (type_of_phases[i] == D2_HISTOGRAM) {
	cout << "\t" << i << "-th phase is of D2_HISTORAM" << endl;
      } else if (type_of_phases[i] == D2_EUCLIDEAN_L2) {
	cout << "\t" << i << "-th phase is of D2_Euclidean_L2 " << endl;
      } else if (type_of_phases[i] == D2_WORD_EMBED) {
	cout << "\t" << i << "-th phase is of D2_WORD_EMBED" << endl;
      }
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
			&dimension_of_phases[0],
			&type_of_phases[0]);

  if (num_of_batches == 0 && is_pre_processed) num_of_batches = 1;

  if (err == 0 && num_of_batches == 0) {  
    d2_read(filename, meta_filename, &data);  
  } else if (num_of_batches > 0 && world_rank == 0) {
    d2_read(filename, meta_filename, &data);      
    d2_write_split(filename, &data, num_of_batches, is_pre_processed);    
    d2_free(&data);
    if (world_rank == 0) {  cout << "[Finish!]" <<endl; }
#ifdef __USE_MPI__
  MPI_Finalize();
#endif
    return 0;
  } else if (err != 0) {
    cerr << "Allocation Failed!" << endl;
  }
  
  /* data structure storing information about centroids of clusters */
  mph c; 
  c.ph = NULL; // make sure c is (re-)initialized 
  std::string name_hashValue;

  /* clustering data and obtain the centroids */
  if (!is_eval) { 
    if (world_rank == 0) {
      if (selected_phase >= 0 && size_of_phases > 1) {
	cout << "Clustering upon " << selected_phase <<"-th phase" << endl;
      } else if (selected_phase < 0 && size_of_phases > 1) {
	cout << "Clustering upon all phases (more than one)" << endl;
      }
    }


    srand (time(NULL));
    int hashNumber=rand() % 1000000;
#ifdef __USE_MPI__
    MPI_Bcast(&hashNumber, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    std::stringstream ss; ss << hashNumber;
    name_hashValue = std::string(std::string(filename) + "_" + ss.str());


    if (is_load) {
      assert(centroid_filename);
      data.num_of_labels = number_of_clusters;
      d2_init_centroid(&data, &c, selected_phase, true);
      d2_read(centroid_filename, meta_filename, &c);        
    }

    if (number_of_clusters == 1) {
      if (d2_alg_type == D2_CENTROID_BADMM) {
	extern BADMM_options *p_badmm_options;
	extern BADMM_options badmm_cen_options;
	p_badmm_options = &badmm_cen_options;
      }	else if (d2_alg_type == D2_CENTROID_GRADDEC) {
	// TBA
      } else if (d2_alg_type == D2_CENTROID_ADMM) {
	// TBA
      }
    }

    d2_clustering(number_of_clusters, 
		  max_iters, 
		  &data, 
		  &c, 
		  selected_phase,
		  use_triangle,
		  name_hashValue.c_str());

    if (world_rank == 0) {
      if (output_filename) name_hashValue = std::string(output_filename);
      d2_write((name_hashValue + "_c.d2").c_str(), &c);
    }
  } 
  
  if (is_eval) {
    assert(centroid_filename);
    std::string cf=std::string(centroid_filename);
    d2_assignment(number_of_clusters, 
		  &data, 
		  &c, 
		  selected_phase, 
		  centroid_filename,
		  meta_filename);
    cf=cf.substr(0, cf.rfind("_c")); // centroid filename should be in match with *_c.d2?
    size_t rpos = cf.rfind("_");
    if (rpos == std::string::npos) {
      // if self provided centroid file is used
      name_hashValue = std::string(filename);
    } else {
      // retrieve the randomized running ID
      cf=cf.substr(rpos);
      name_hashValue = std::string(filename)+ cf;
    }
  }

  if (output_filename) name_hashValue = std::string(output_filename);
  d2_write_labels((name_hashValue + ".label").c_str(), &data);
  d2_write_labels_serial((std::string(filename) + ".ind").c_str(), 
			 name_hashValue.c_str(), &data);

  d2_free(&data);
  d2_free(&c);

  if (world_rank == 0) {  cout << "[Finish!]" <<endl; }
#ifdef __USE_MPI__
  MPI_Finalize();
#endif
  return 0;
}

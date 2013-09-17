#include <fstream>
#include <iostream>
#include <vector>

#include "d2_clustering.h"

#define MODAL_SIZE (2)
using std::vector;
using std::ifstream;
using std::cout;
using std::endl;


int main(int argc, char *argv[])
{ 

  int size_of_modalities = MODAL_SIZE;
  int dimensions_of_modalities[MODAL_SIZE];

  vector<int> size_of_supports;
  vector<SCALAR> data_block_supp, data_block_w;

  ifstream fp("mountaindat.txt");
  int istart = 0;
  while (fp.good()) {
    for (int i=0; i<size_of_modalities; ++i) {
      int dim, size;
      fp >> dim >> size;
      if (0 == istart) dimensions_of_modalities[i] = dim;
      size_of_supports.push_back(size);

      SCALAR c;
      for (int j=0; j<size; ++j) {fp >> c; data_block_w.push_back(c);}
      for (int j=0; j<size * dim; ++j) {fp >> c; data_block_supp.push_back(c);}
      
    }
    ++istart;
  }
  fp.close();

  int size_of_samples = istart;

  //////////////////////////////////////
  d2_initialize(size_of_modalities, dimensions_of_modalities);
  d2_assign_data(size_of_samples, &size_of_supports[0], &data_block_supp[0], &data_block_w[0]);

  return 0;
}

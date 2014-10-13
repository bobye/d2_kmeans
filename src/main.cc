#include <stdio.h>
#include <string>
#include <vector>
#include <iostream>


#include "d2_clustering.h"
#include "util.hh"


// ./bin/test mountaindat.txt 2 1000 3,3 6,11
int main(int argc, char *argv[])
{ 
  using namespace std;

  int size_of_phases = atoi(argv[2]);
  int size_of_samples = atoi(argv[3]);

  vector<int> dimension_of_phases(size_of_phases);  
  vector<int> avg_strides(size_of_phases);

  vector<string> ss1 = split(string(argv[4]), ',');
  vector<string> ss2 = split(string(argv[5]), ',');

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
    fp = fopen(argv[1], "r+");
    d2_load(fp, &data);  
    fclose(fp);
  } else {
    cerr << "Allocation Failed!" << endl;
  }

  sph c;
  d2_centroid_sphBregman(&data, 0, NULL, &c);
  

  d2_free(&data);

  return 0;
}

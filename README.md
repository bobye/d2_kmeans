d2_clustering
=============

scalable and parallel D2 (discrete distribution) clustering

- Algorithm Prototypes
  matlab/
    d2_clustering.m --- old scripts
    d2clusters.m --- old clustering framework
    profile_centroid.m --- profiling the convergence of centroid updates
    profile_kantorovich.m --- profiling LP solution of transportation problem
    centroid_sph*.m --- computing centroid of a single phase

- C/C++ sources (algorithms should be implemented in C98, but tests or others can be in C++)
  include/
  src/
  Makefile

- Dataset

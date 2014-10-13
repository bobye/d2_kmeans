#ifndef _D2_CLUSTERING_H_
#define _D2_CLUSTERING_H_

#ifdef __cplusplus
extern "C" {
#endif

#define SCALAR double
#define SCALAR_SCANF_TYPE ("%lf")

  typedef struct {
    int avg_str;
    int dim;
    int *p_str;
    SCALAR *p_supp;
    SCALAR *p_w;  
  } sph;
  

  typedef struct {
    int s_ph; // size of phases
    int size;
    sph *ph;
  } mph;

  int d2_allocate(mph *p_data,
		  const int size_of_phases,
		  const int size_of_samples,
		  const int *avg_strides,
		  const int *dimension_of_phases);

  int d2_load(void *fp, mph *p_data);

  int d2_free(mph *p_data);

#ifdef __cplusplus
}
#endif

#endif /* _D2_CLUSTERING_H_ */

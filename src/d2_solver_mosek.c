#include "d2_solver.h"
#include "mosek.h"
#include <stdio.h>

#define SCALAR double

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_FREE(x)         free(x)

#elif defined __USE_MKL__
#include <mkl.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) mkl_malloc( (x) *sizeof(SCALAR), 16) 
#define _D2_MALLOC_INT(x)       (int *) mkl_malloc( (x) *sizeof(int), 16)
#define _D2_FREE(x)         mkl_free(x)

#elif defined __GNUC__
#include <cblas.h>
#include <lapacke.h>
#define _D2_MALLOC_SCALAR(x)       (SCALAR *) malloc( (x) *sizeof(SCALAR)) 
#define _D2_MALLOC_INT(x)       (int *) malloc( (x) *sizeof(int))
#define _D2_FREE(x)         free(x)

#endif

MSKenv_t   env = NULL;
/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void *handle,
                            MSKCONST char str[])
{
  printf("%s",str);
} /* printstr */


void d2_solver_setup() {
  MSKrescodee r;

  /* Create the mosek environment. */
  r = MSK_makeenv(&env, NULL);

  if ( r != MSK_RES_OK ) {
      /* In case of an error print error code and description. */      
      char symname[MSK_MAX_STR_LEN];
      char desc[MSK_MAX_STR_LEN];
      
      printf("An error occurred while optimizing.\n");     
      MSK_getcodedesc (r,
                       symname,
                       desc);
      printf("Error %s - '%s'\n",symname,desc);
  }
  return;
}

void d2_solver_debug() {
  int i,j;
  double X[] = {83.718147, 0.355520, 2.771609, 47.223366, 32.704613, 24.236309, 44.592049, -23.987589, -70.987701, 51.051491, -0.224495, -6.329401, 83.664337, 10.054308, 25.707617, 66.965271, -31.481861, -46.937881, 71.953934, -12.366786, -18.309299, 69.359741, 16.177128, 25.328596};
  double Y[] = {80.360435, -12.656796, -23.847801, 14.429349, 10.828419, -4.769543, 27.646425, -7.103071, -2.075422, 52.567692, -27.300863, -27.065271};
  double wX[] = {0.210949, 0.017211, 0.045296, 0.085237, 0.192448, 0.151259, 0.123296, 0.174304};
  double wY[] = {0.546734, 0.062044, 0.076032, 0.315190};

  double x[32], lambda[12];

  d2_match_by_coordinates(3, 8, X, wX, 4, Y, wY, x, lambda);
  //for (i=0; i<12; ++i) printf("%f ", lambda[i]);
}

void d2_solver_release() {
  MSK_deleteenv(&env);
}

double d2_match_by_sample(int d, int n, double *X, int m, double *Y, 
			  /** OUT **/ double *x, /** OUT **/ double *lambda) {
  int i;
  double *wX, *wY, fval;
  wX = _D2_MALLOC_SCALAR(n);
  wY = _D2_MALLOC_SCALAR(m);
  for (i=0; i<n; ++i) wX[i] = 1./n;
  for (i=0; i<m; ++i) wY[i] = 1./m;

  fval = d2_match_by_coordinates(d, n, X, wX, m, Y, wY, x, lambda);
  _D2_FREE(wX);
  _D2_FREE(wY);

  return fval;
}

double d2_match_by_distmat(int n, int m, double *C, double *wX, double *wY,
			   /** OUT **/ double *x, /** OUT **/ double *lambda) {
  
  const MSKint32t numvar = n * m,
                  numcon = n + m;
  MSKtask_t    task = NULL;
  MSKrescodee r;
  MSKint32t    i,j;
  double ones[] = {1.0, 1.0}, fval = 0.0;
  MSKint32t *asub = (MSKint32t *) malloc(2*m*n* sizeof(MSKint32t));

  for (j=0; j<m; ++j) {
    for (i=0; i<n; ++i) {
      asub[2*(i +j*n)] = i;
      asub[2*(i +j*n) +1] = n+j;
    }
  }

  /* Create the optimization task. */
  r = MSK_maketask(env,numcon,numvar,&task);
  if ( r==MSK_RES_OK ) {
    //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
  }

  r = MSK_appendcons(task,numcon);
  r = MSK_appendvars(task,numvar); 

  for (j=0; j<numvar && r == MSK_RES_OK; ++j) {
    r = MSK_putcj(task,j,C[j]);
    r = MSK_putvarbound(task, 
			j,           /* Index of variable.*/
			MSK_BK_LO,      /* Bound key.*/
			0.0,      /* Numerical value of lower bound.*/
                        +MSK_INFINITY);     /* Numerical value of upper bound.*/

    r = MSK_putacol(task, 
		    j,           /* Index of variable.*/
		    2,           /* Number of non-zeros in column j.*/
		    asub+j*2,
		    ones);
  }

  for (i=0; i<n && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(task,
			i,
			MSK_BK_FX,
			wX[i],
			wX[i]);
  for (i=0; i<m && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(task,
			i+n,
			MSK_BK_FX,
			wY[i],
			wY[i]);

  if (r == MSK_RES_OK) {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }

  if ( r==MSK_RES_OK )
    {
      MSKrescodee trmcode;
      /* Run optimizer */
      r = MSK_optimizetrm(task,&trmcode);      

      /* Print a summary containing information
         about the solution for debugging purposes. */
      //MSK_solutionsummary (task,MSK_STREAM_LOG);      
      if ( r==MSK_RES_OK ) {
	MSKsolstae solsta;
	if ( r==MSK_RES_OK )
          r = MSK_getsolsta (task,
                             MSK_SOL_BAS,
                             &solsta);
	switch(solsta)
	{
          case MSK_SOL_STA_OPTIMAL:   
          case MSK_SOL_STA_NEAR_OPTIMAL:
          {
            //double *xx = (double*) calloc(numvar,sizeof(double));
	    MSK_getprimalobj(task, MSK_SOL_BAS, &fval);
            if ( x )
            {
              MSK_getxx(task,
                        MSK_SOL_BAS,    /* Request the basic solution. */
                        x);        
              //printf("Optimal primal solution\n");
              //for(j=0; j<numvar; ++j) printf("x[%d]: %e\n",j,xx[j]);

              //free(x);
            }

	    if (lambda) 
	    {
	      MSK_gety (task,
                        MSK_SOL_BAS,    /* Request the dual solution: be careful about exact +- of variables */
                        lambda);
	    }	    
            else 
              r = MSK_RES_ERR_SPACE;

            break;
          }
          case MSK_SOL_STA_DUAL_INFEAS_CER:
          case MSK_SOL_STA_PRIM_INFEAS_CER:
          case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
          case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:  
            printf("Primal or dual infeasibility certificate found.\n");
            break;
          case MSK_SOL_STA_UNKNOWN:
          {
            char symname[MSK_MAX_STR_LEN];
            char desc[MSK_MAX_STR_LEN];

            /* If the solutions status is unknown, print the termination code
               indicating why the optimizer terminated prematurely. */
            
            MSK_getcodedesc(trmcode,
                            symname,
                            desc);
            
            printf("The solution status is unknown.\n");
            printf("The optimizer terminitated with code: %s\n",symname);
            break;
          }
          default:
            printf("Other solution status.\n");
            break;
        }
      }
    }
  return fval;
}


double d2_match_by_coordinates(int d, int n, double *X, double *wX, int m, double *Y,  double*wY,
			       /** OUT **/ double *x, /** OUT **/ double *lambda) {
  double *C, fval;
  C = _D2_MALLOC_SCALAR(n*m);
  _dpdist2(d, n, m, X, Y, C);
  fval = d2_match_by_distmat(n, m, C, wX, wY, x, lambda);
  _D2_FREE(C);
  return fval;
}

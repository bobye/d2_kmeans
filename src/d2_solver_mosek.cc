#include "d2_solver.h"
#include "mosek.h"
#include <stdio.h>
#include "blas_util.h"

static MSKenv_t   env = NULL;
//static MSKtask_t *task_seq = NULL;
//static size_t task_seq_size = 0;

#include <utility>  
#include <map>
using std::pair;
using std::make_pair;
using std::map;
static map< pair<int, int>, MSKtask_t > task_mapper;

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
  /*
  task_seq = (MSKtask_t*) malloc (num * sizeof(MSKtask_t));
  for (i=0; i<num; ++i) task_seq[i] = NULL;
  task_seq_size = num;
  */
  return;
}

void d2_solver_release() {
  /*
  size_t i;
  for (i=0; i<task_seq_size; ++i) 
    if (task_seq[i] != NULL) MSK_deletetask(&task_seq [i]);
  */
  for (map< pair<int, int>, MSKtask_t >::iterator it=task_mapper.begin(); it!=task_mapper.end(); ++it) MSK_deletetask(&(it->second));
  MSK_deleteenv(&env);
}


double d2_match_by_distmat(int n, int m, SCALAR *C, SCALAR *wX, SCALAR *wY,
			   __OUT__ SCALAR *x, __OUT__ SCALAR *lambda, size_t index) {
  
  const MSKint32t numvar = n * m,
                  numcon = n + m - 1;
  MSKtask_t    *p_task;
  MSKrescodee r = MSK_RES_OK;
  MSKint32t    i,j;
  double fval = 0.0;

  /* check if it is in the mode of multiple phase or single phase */
  //  p_task = &task_seq[index];

  if (task_mapper.find(make_pair (n, m)) == task_mapper.end()) {
    task_mapper[make_pair (n, m)] = NULL;
  }
  p_task = &task_mapper[make_pair (n, m)];

  if (*p_task == NULL) {
  MSKint32t *asub;
  double ones[2] = {1.0, 1.0};
  asub = (MSKint32t *) malloc(2*m*n* sizeof(MSKint32t));
  for (j=0; j<m; ++j) {
    for (i=0; i<n; ++i) {
      asub[2*(i +j*n)] = i;
      asub[2*(i +j*n) +1] = n+j;
    }
  }

  /* Create the optimization task. */
  r = MSK_maketask(env,numcon,numvar,p_task); 
  if ( r==MSK_RES_OK ) {
    //r = MSK_linkfunctotaskstream(*p_task,MSK_STREAM_LOG,NULL,printstr);
  }

  r = MSK_appendcons(*p_task,numcon);
  r = MSK_appendvars(*p_task,numvar); 

  for (j=0; j<numvar && r == MSK_RES_OK; ++j) {
    r = MSK_putvarbound(*p_task, 
			j,           /* Index of variable.*/
			MSK_BK_LO,      /* Bound key.*/
			0.0,      /* Numerical value of lower bound.*/
                        +MSK_INFINITY);     /* Numerical value of upper bound.*/

    i = (j >= (m-1)*n) ? 1 : 2;
    r = MSK_putacol(*p_task, 
		    j,           /* Index of variable.*/
		    i,           /* Number of non-zeros in column j.*/
		    asub+j*2,
		    ones);
  }
  free(asub);

  if (r == MSK_RES_OK) {
    r = MSK_putobjsense(*p_task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }
  /* Disable presolve: may lead to minor improvement */
  // r = MSK_putintparam(task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF);
  /* set network flow problem */
  r = MSK_putintparam(*p_task, MSK_IPAR_OPTIMIZER,  MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX);
  /* disable multi-threads */
  r = MSK_putintparam(*p_task, MSK_IPAR_NUM_THREADS, 1);
  } else {
  }

  /* modify an existing task and re-optimize */
  for (j=0; j<numvar && r == MSK_RES_OK; ++j) {
    r = MSK_putcj(*p_task,j,C[j]);
  }

  for (i=0; i<n && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(*p_task,
			i,
			MSK_BK_FX,
			wX[i],
			wX[i]);
  for (i=0; i<m-1 && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(*p_task,
			i+n,
			MSK_BK_FX,
			wY[i],
			wY[i]);


  if ( r==MSK_RES_OK )
    {
      MSKrescodee trmcode;

      /* Run optimizer */
      r = MSK_optimizetrm(*p_task,&trmcode);      

      /* Print a summary containing information
         about the solution for debugging purposes. */
      //MSK_solutionsummary (task,MSK_STREAM_LOG);      
      if ( r==MSK_RES_OK ) {
	MSKsolstae solsta;
	r = MSK_getsolsta (*p_task,
			   MSK_SOL_BAS,
			   &solsta);
	switch(solsta)
	{
          case MSK_SOL_STA_OPTIMAL:   
          case MSK_SOL_STA_NEAR_OPTIMAL:
          {
            //double *xx = (double*) calloc(numvar,sizeof(double));
	    MSK_getprimalobj(*p_task, MSK_SOL_BAS, &fval);
            if ( x )
            {
#ifdef _D2_DOUBLE
              MSK_getxx(*p_task,
                        MSK_SOL_BAS,    /* Request the basic solution. */
                        x);        
#endif
              //printf("Optimal primal solution\n");
              //for(j=0; j<numvar; ++j) printf("x[%d]: %e\n",j,xx[j]);

              //free(x);
            }

	    if (lambda) 
	    {
#ifdef _D2_DOUBLE
	      MSK_gety (*p_task,
                        MSK_SOL_BAS,    /* Request the dual solution: be careful about exact +- of variables */
                        lambda);
#endif
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

  //  MSK_deletetask(&task);

  return fval;
}




/**
 * Main codes ends and extra codes begins.
 */

double d2_match_by_distmat_qp(int n, int m, 
			      SCALAR *C, SCALAR *L, SCALAR rho, 
			      SCALAR *lw, SCALAR *rw, 
			      SCALAR *x0, 
			      /** OUT **/ SCALAR *x) {
  const MSKint32t numvar = n*m + n, numcon = n + m;
  MSKtask_t task = NULL;
  MSKrescodee r;
  MSKint32t i, j;
  double fval = 0.0, ones[] = {1.0, 1.0};
  MSKint32t *asub = (MSKint32t *) malloc(2*m*n*sizeof(MSKint32t));
  MSKint32t *bsub = (MSKint32t *) malloc(n*sizeof(MSKint32t));
  MSKint32t *qsub = (MSKint32t *) malloc(n*sizeof(MSKint32t));
  MSKrealt *qval = (MSKrealt *) malloc(n*sizeof(MSKrealt));
  MSKrealt *xx = (MSKrealt *) malloc(numvar*sizeof(MSKrealt));
  for (j=0; j<m; ++j)
    for (i=0; i<n; ++i) {
      asub[2*(i+j*n)] = i;
      asub[2*(i+j*n) + 1] = n+j;
    }
  for (j=0; j<n; ++j) {
    bsub[j] = j;
    qsub[j] = n*m + j;
    qval[j] = rho;
  }
  /* Create the optimization task. */
  r = MSK_maketask(env,numcon,numvar,&task);
  if ( r==MSK_RES_OK ) {
    //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
  }

  r = MSK_appendcons(task,numcon);
  r = MSK_appendvars(task,numvar); 

  for (j=0; j<n*m && r == MSK_RES_OK; ++j) {
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
  for (j=0; j<n && r == MSK_RES_OK; ++j) {
    r = MSK_putvarbound(task,
			n*m + j,
			MSK_BK_FR,
			-MSK_INFINITY,
			+MSK_INFINITY);
    r = MSK_putacol(task,
		    n*m + j,
		    1,
		    bsub + j,
		    ones);
  }

  for (i=0; i<n && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(task,
			i,
			MSK_BK_FX,
			lw[i] - L[i],
			lw[i] - L[i]);
  for (i=0; i<m && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(task,
			i+n,
			MSK_BK_FX,
			rw[i],
			rw[i]);

  if (r == MSK_RES_OK) {
    r = MSK_putqobj(task,n,qsub,qsub,qval); 
  }


  if (r == MSK_RES_OK) {
    r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  }

  if (r == MSK_RES_OK) 
    {
      MSKrescodee trmcode;
      r = MSK_optimizetrm(task, &trmcode);

      if (r == MSK_RES_OK) {
	MSKsolstae solsta;
	r = MSK_getsolsta( task,
			   MSK_SOL_ITR,  /* Request the interior solution. */ 
			   &solsta);
	switch(solsta)
	  {
	  case MSK_SOL_STA_OPTIMAL:
	  case MSK_SOL_STA_NEAR_OPTIMAL:
	    MSK_getprimalobj(task, MSK_SOL_ITR, &fval);
	    if ( x ) {
	      MSK_getxx(task,
			MSK_SOL_ITR,
			xx);
	      for (i=0; i<n*m; ++i) x[i] = xx[i];	     
	    }
	    break;
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
  
  MSK_deletetask(&task); 
  free(asub); free(bsub); free(qsub); free(qval); free(xx);
  return fval;
}


double d2_qpsimple(int n, int count, SCALAR *c, /** OUT **/ SCALAR *w) {
  const MSKint32t numvar = n, numcon = 1;
  MSKtask_t task = NULL;
  MSKrescodee r;
  MSKint32t j;
  MSKint32t zeroInt = 0;
  double one = 1, fval = 0.0;

  MSKint32t *bsub = (MSKint32t *) malloc(n*sizeof(MSKint32t));
  MSKint32t *qsub = (MSKint32t *) malloc(n*sizeof(MSKint32t));
  MSKrealt *qval = (MSKrealt *) malloc(n*sizeof(MSKrealt));

  for (j=0; j<n; ++j) {
    qsub[j] = j; // Q index
    qval[j] = count; // Q value
  }
  r=MSK_maketask(env, numcon, numvar, &task);
  if ( r==MSK_RES_OK ) {
    //r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
  }

  r=MSK_appendcons(task, numcon);
  r=MSK_appendvars(task, numvar);
  for (j=0; j<n; ++j) {
    r = MSK_putcj(task, j, -c[j]);
    r = MSK_putvarbound(task, j, MSK_BK_LO, 0.0, +MSK_INFINITY);
    r = MSK_putacol(task, j, 1, &zeroInt, &one);
  }
  r = MSK_putconbound(task, 0, MSK_BK_FX, 1.0, 1.0);
  r = MSK_putqobj(task, n, qsub, qsub, qval);
  r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);


  if (r == MSK_RES_OK) {
      MSKrescodee trmcode;
      r = MSK_optimizetrm(task, &trmcode);
      if (r == MSK_RES_OK) {
	MSKsolstae solsta;
	r = MSK_getsolsta( task,
			   MSK_SOL_ITR,  /* Request the interior solution. */ 
			   &solsta);
	switch(solsta)
	  {
	  case MSK_SOL_STA_OPTIMAL:
	  case MSK_SOL_STA_NEAR_OPTIMAL:
	    MSK_getprimalobj(task, MSK_SOL_ITR, &fval);
	    if ( w ) {
#ifdef _D2_DOUBLE
	      MSK_getxx(task,
			MSK_SOL_ITR,
			w);
#endif
	    }
	    break;
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

  MSK_deletetask(&task); 
  free(bsub); free(qsub); free(qval); 

  return fval;
}

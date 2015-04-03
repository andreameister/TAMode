/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2010/12/01 22:21:04 $
 * -----------------------------------------------------------------
 * Programmer(s): Alan C. Hindmarsh and Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the optional input and output
 * functions for the CVODE solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"

#define ZERO (0.0)
#define ONE  (1.0)

#define lrw (cv_mem->cv_lrw)
#define liw (cv_mem->cv_liw)
#define lrw1 (cv_mem->cv_lrw1)
#define liw1 (cv_mem->cv_liw1)
#define TRUE 1
#define FALSE 0

/* 
 * =================================================================
 * CVODE optional input functions
 * =================================================================
 */

/* 
 * CVodeSetErrHandlerFn
 *
 * Specifies the error handler function
 */

int CVodeSetErrHandlerFn(void * const cvode_mem, CVErrHandlerFn ehfun, void * const eh_data)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetErrHandlerFn", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_ehfun = ehfun;
  cv_mem->cv_eh_data = eh_data;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetErrFile
 *
 * Specifies the FILE pointer for output (NULL means no messages)
 */

int CVodeSetErrFile(void * const cvode_mem, FILE * const errfp)
{
  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetErrFile", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  CVodeMem cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_errfp = errfp;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetUserData
 *
 * Specifies the user data pointer for f
 */

int CVodeSetUserData(void * const cvode_mem, void * const user_data)
{
    if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetUserData", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  CVodeMem cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_user_data = user_data;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxNumSteps
 *
 * Specifies the maximum number of integration steps
 */

int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetMaxNumSteps", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  /* Passing mxsteps=0 sets the default. Passing mxsteps<0 disables the test. */

  if (mxsteps == 0)
    cv_mem->cv_mxstep = MXSTEP_DEFAULT;
  else
    cv_mem->cv_mxstep = mxsteps;

  return(CV_SUCCESS);
}

/* 
 * CVodeSetMaxHnilWarns
 *
 * Specifies the maximum number of warnings for small h
 */

int CVodeSetMaxHnilWarns(void *cvode_mem, int mxhnil)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeSetMaxHnilWarns", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  cv_mem->cv_mxhnil = mxhnil;

  return(CV_SUCCESS);
}




/* 
 * =================================================================
 * CVODE optional output functions
 * =================================================================
 */

/* 
 * Readability constants
 */

#define nst            (cv_mem->cv_nst)
#define nfe            (cv_mem->cv_nfe)
#define ncfn           (cv_mem->cv_ncfn)
#define netf           (cv_mem->cv_netf)
#define nni            (cv_mem->cv_nni)
#define nsetups        (cv_mem->cv_nsetups)
#define qu             (cv_mem->cv_qu)
#define next_q         (cv_mem->cv_next_q)
#define ewt            (cv_mem->cv_ewt)  
#define hu             (cv_mem->cv_hu)
#define next_h         (cv_mem->cv_next_h)
#define h0u            (cv_mem->cv_h0u)
#define tolsf          (cv_mem->cv_tolsf)  
#define acor           (cv_mem->cv_acor)
#define lrw            (cv_mem->cv_lrw)
#define liw            (cv_mem->cv_liw)
#define nge            (cv_mem->cv_nge)
#define iroots         (cv_mem->cv_iroots)
#define nor            (cv_mem->cv_nor)
#define sldeton        (cv_mem->cv_sldeton)
#define tn             (cv_mem->cv_tn)
#define efun           (cv_mem->cv_efun)


/* 
 * CVodeGetIntegratorStats
 *
 * Returns integrator statistics
 */

int CVodeGetIntegratorStats(void *cvode_mem, long int *nsteps, long int *nfevals, 
                            long int *nlinsetups, long int *netfails, int *qlast, 
                            int *qcur, double *hinused, double *hlast, 
                            double *hcur, double *tcur)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeGetIntegratorStats", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nsteps = nst;
  *nfevals = nfe;
  *nlinsetups = nsetups;
  *netfails = netf;
  *qlast = qu;
  *qcur = next_q;
  *hinused = h0u;
  *hlast = hu;
  *hcur = next_h;
  *tcur = tn;

  return(CV_SUCCESS);
}


/* 
 * CVodeGetNonlinSolvStats
 *
 * Returns nonlinear solver statistics
 */

int CVodeGetNonlinSolvStats(void *cvode_mem, long int *nniters, 
                            long int *nncfails)
{
  CVodeMem cv_mem;

  if (cvode_mem==NULL) {
    CVProcessError(NULL, CV_MEM_NULL, "CVODE", "CVodeGetNonlinSolvStats", MSGCV_NO_MEM);
    return(CV_MEM_NULL);
  }

  cv_mem = (CVodeMem) cvode_mem;

  *nniters = nni;
  *nncfails = ncfn;

  return(CV_SUCCESS);
}

/*-----------------------------------------------------------------*/

char *CVodeGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(24*sizeof(char));

  switch(flag) {
  case CV_SUCCESS:
    sprintf(name,"CV_SUCCESS");
    break;
  case CV_TSTOP_RETURN:
    sprintf(name,"CV_TSTOP_RETURN");
    break;
  case CV_ROOT_RETURN:
    sprintf(name,"CV_ROOT_RETURN");
    break;
  case CV_TOO_MUCH_WORK:
    sprintf(name,"CV_TOO_MUCH_WORK");
    break;
  case CV_TOO_MUCH_ACC:
    sprintf(name,"CV_TOO_MUCH_ACC");
    break;
  case CV_ERR_FAILURE:
    sprintf(name,"CV_ERR_FAILURE");
    break;
  case CV_CONV_FAILURE:
    sprintf(name,"CV_CONV_FAILURE");
    break;
  case CV_LINIT_FAIL:
    sprintf(name,"CV_LINIT_FAIL");
    break;
  case CV_LSETUP_FAIL:
    sprintf(name,"CV_LSETUP_FAIL");
    break;
  case CV_LSOLVE_FAIL:
    sprintf(name,"CV_LSOLVE_FAIL");
    break;
  case CV_RHSFUNC_FAIL:
    sprintf(name,"CV_RHSFUNC_FAIL");
    break;
  case CV_FIRST_RHSFUNC_ERR:
    sprintf(name,"CV_FIRST_RHSFUNC_ERR");
    break;
  case CV_REPTD_RHSFUNC_ERR:
    sprintf(name,"CV_REPTD_RHSFUNC_ERR");
    break;
  case CV_UNREC_RHSFUNC_ERR:
    sprintf(name,"CV_UNREC_RHSFUNC_ERR");
    break;
  case CV_RTFUNC_FAIL:
    sprintf(name,"CV_RTFUNC_FAIL");
    break;
  case CV_MEM_FAIL:
    sprintf(name,"CV_MEM_FAIL");
    break;
  case CV_MEM_NULL:
    sprintf(name,"CV_MEM_NULL");
    break;
  case CV_ILL_INPUT:
    sprintf(name,"CV_ILL_INPUT");
    break;
  case CV_NO_MALLOC:
    sprintf(name,"CV_NO_MALLOC");
    break;
  case CV_BAD_K:
    sprintf(name,"CV_BAD_K");
    break;
  case CV_BAD_T:
    sprintf(name,"CV_BAD_T");
    break;
  case CV_BAD_DKY:
    sprintf(name,"CV_BAD_DKY");
    break;
  case CV_TOO_CLOSE:
    sprintf(name,"CV_TOO_CLOSE");
    break;    
  default:
    sprintf(name,"NONE");
  }

  return(name);
}


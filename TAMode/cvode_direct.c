/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2010/12/01 22:21:04 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for the CVDLS linear solvers
 * -----------------------------------------------------------------
 */

/* 
 * =================================================================
 * IMPORTED HEADER FILES
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_impl.h"
#include "cvode_direct_impl.h"

/* 
 * =================================================================
 * FUNCTION SPECIFIC CONSTANTS
 * =================================================================
 */

/* Constant for DQ Jacobian approximation */
#define MIN_INC_MULT (1000.0)

#define ZERO         (0.0)
#define ONE          (1.0)
#define TWO          (2.0)

/*
 * =================================================================
 * READIBILITY REPLACEMENTS
 * =================================================================
 */

#define f         (cv_mem->cv_f)
#define user_data (cv_mem->cv_user_data)
#define uround    (cv_mem->cv_uround)
#define nst       (cv_mem->cv_nst)
#define tn        (cv_mem->cv_tn)
#define h         (cv_mem->cv_h)
#define gamma     (cv_mem->cv_gamma)
#define gammap    (cv_mem->cv_gammap)
#define gamrat    (cv_mem->cv_gamrat)
#define ewt       (cv_mem->cv_ewt)

#define lmem      (cv_mem->cv_lmem)

#define mtype     (cvdls_mem->d_type)
#define n         (cvdls_mem->d_n)
#define ml        (cvdls_mem->d_ml)
#define mu        (cvdls_mem->d_mu)
#define smu       (cvdls_mem->d_smu)
#define jacDQ     (cvdls_mem->d_jacDQ)
#define djac      (cvdls_mem->d_djac)
#define bjac      (cvdls_mem->d_bjac)
#define M         (cvdls_mem->d_M)
#define nje       (cvdls_mem->d_nje)
#define nfeDQ     (cvdls_mem->d_nfeDQ)
#define last_flag (cvdls_mem->d_last_flag)

/* 
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */
              
/*
 * CVDlsSetDenseJacFn specifies the dense Jacobian function.
 */
int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsSetDenseJacFn", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsSetDenseJacFn", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  if (jac != NULL) {
    jacDQ = FALSE;
    djac = jac;
  } else {
    jacDQ = TRUE;
  }

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetNumJacEvals returns the number of Jacobian evaluations.
 */
int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetNumJacEvals", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetNumJacEvals", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *njevals = nje;

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetNumRhsEvals returns the number of calls to the ODE function
 * needed for the DQ Jacobian approximation.
 */
int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetNumRhsEvals", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetNumRhsEvals", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *nfevalsLS = nfeDQ;

  return(CVDLS_SUCCESS);
}

/*
 * CVDlsGetReturnFlagName returns the name associated with a CVDLS
 * return value.
 */
char *CVDlsGetReturnFlagName(long int flag)
{
  char *name;

  name = (char *)malloc(30*sizeof(char));

  switch(flag) {
  case CVDLS_SUCCESS:
    sprintf(name,"CVDLS_SUCCESS");
    break;   
  case CVDLS_MEM_NULL:
    sprintf(name,"CVDLS_MEM_NULL");
    break;
  case CVDLS_LMEM_NULL:
    sprintf(name,"CVDLS_LMEM_NULL");
    break;
  case CVDLS_ILL_INPUT:
    sprintf(name,"CVDLS_ILL_INPUT");
    break;
  case CVDLS_MEM_FAIL:
    sprintf(name,"CVDLS_MEM_FAIL");
    break;
  case CVDLS_JACFUNC_UNRECVR:
    sprintf(name,"CVDLS_JACFUNC_UNRECVR");
    break;
  case CVDLS_JACFUNC_RECVR:
    sprintf(name,"CVDLS_JACFUNC_RECVR");
    break;
  default:
    sprintf(name,"NONE");
  }

  return(name);
}

/*
 * CVDlsGetLastFlag returns the last flag set in a CVDLS function.
 */
int CVDlsGetLastFlag(void *cvode_mem, long int *flag)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDLS_MEM_NULL, "CVDLS", "CVDlsGetLastFlag", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  if (lmem == NULL) {
    CVProcessError(cv_mem, CVDLS_LMEM_NULL, "CVDLS", "CVDlsGetLastFlag", MSGD_LMEM_NULL);
    return(CVDLS_LMEM_NULL);
  }
  cvdls_mem = (CVDlsMem) lmem;

  *flag = last_flag;

  return(CVDLS_SUCCESS);
}

/* 
 * =================================================================
 * DQ JACOBIAN APPROXIMATIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * cvDlsDenseDQJac 
 * -----------------------------------------------------------------
 * This routine generates a dense difference quotient approximation to
 * the Jacobian of f(t,y). It assumes that a dense matrix of type
 * DlsMat is stored column-wise, and that elements within each column
 * are contiguous. The address of the jth column of J is obtained via
 * the macro DENSE_COL and this pointer is associated with an N_Vector
 * using the N_VGetArrayPointer/N_VSetArrayPointer functions. 
 * Finally, the actual computation of the jth column of the Jacobian is 
 * done with a call to N_VLinearSum.
 * -----------------------------------------------------------------
 */ 

int cvDlsDenseDQJac(long int N, double t,
                    N_Vector y, N_Vector fy, 
                    DlsMat Jac, void *data,
                    N_Vector tmp1, N_Vector tmp2, __attribute__((unused)) N_Vector tmp3)
{
  double fnorm, minInc, inc, inc_inv, yjsaved, srur;
  double *tmp2_data, *y_data, *ewt_data;
  N_Vector ftemp, jthCol;
  long int j;
  int retval = 0;

  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* data points to cvode_mem */
  cv_mem = (CVodeMem) data;
  cvdls_mem = (CVDlsMem) lmem;

  /* Save pointer to the array in tmp2 */
  tmp2_data = N_VGetArrayPointer(tmp2);

  /* Rename work vectors for readibility */
  ftemp = tmp1; 
  jthCol = tmp2;

  /* Obtain pointers to the data for ewt, y */
  ewt_data = N_VGetArrayPointer(ewt);
  y_data   = N_VGetArrayPointer(y);

  /* Set minimum increment based on uround and norm of f */
  srur = RSqrt(uround);
  fnorm = N_VWrmsNorm(fy, ewt);
  minInc = (fnorm != ZERO) ?
           (MIN_INC_MULT * fabs(h) * uround * N * fnorm) : ONE;

  for (j = 0; j < N; j++) {

    /* Generate the jth col of J(tn,y) */

    N_VSetArrayPointer(DENSE_COL(Jac,j), jthCol);

    yjsaved = y_data[j];
    inc = MAX(srur*fabs(yjsaved), minInc/ewt_data[j]);
    y_data[j] += inc;

    retval = f(t, y, ftemp, user_data);
    nfeDQ++;
    if (retval != 0) break;
    
    y_data[j] = yjsaved;

    inc_inv = ONE/inc;
    N_VLinearSum(inc_inv, ftemp, -inc_inv, fy, jthCol);

    DENSE_COL(Jac,j) = N_VGetArrayPointer(jthCol);
  }

  /* Restore original array pointer in tmp2 */
  N_VSetArrayPointer(tmp2_data, tmp2);

  return(retval);
}


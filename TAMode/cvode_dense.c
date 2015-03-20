/*
 * -----------------------------------------------------------------
 * $Revision: 1.12 $
 * $Date: 2010/12/01 22:21:04 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the impleentation file for the CVDENSE linear solver.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "cvode_dense.h"
#include "cvode_direct_impl.h"
#include "cvode_impl.h"

/* Constants */

#define ZERO         (0.0)
#define ONE          (1.0)
#define TWO          (2.0)

/* CVDENSE linit, lsetup, lsolve, and lfree routines */
 
static int cvDenseInit(CVodeMem cv_mem);

static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, int *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector weight,
                        N_Vector ycur, N_Vector fcur);

static void cvDenseFree(CVodeMem cv_mem);



/*
 * -----------------------------------------------------
 * Functions working on DlsMat
 * -----------------------------------------------------
 */
static long int denseGETRF(double **a, unsigned int m, unsigned int n, long int *p) {
    double *col_j, *col_k;
    double temp, mult, a_kj;
    
    /* k-th elimination step number */
    for (unsigned int k=0; k < n; k++) {
        
        col_k  = a[k];
        
        /* find l = pivot row number */
        unsigned int l = k;
        for (unsigned int i=k+1; i < m; i++)
            if (fabs(col_k[i]) > fabs(col_k[l])) l=i;
        p[k] = l;
        
        /* check for zero pivot element */
        if (col_k[l] == ZERO) return(k+1);
        
        /* swap a(k,1:n) and a(l,1:n) if necessary */
        if ( l!= k ) {
            for (unsigned int i=0; i<n; i++) {
                temp = a[i][l];
                a[i][l] = a[i][k];
                a[i][k] = temp;
            }
        }
        
        /* Scale the elements below the diagonal in
         * column k by 1.0/a(k,k). After the above swap
         * a(k,k) holds the pivot element. This scaling
         * stores the pivot row multipliers a(i,k)/a(k,k)
         * in a(i,k), i=k+1, ..., m-1.
         */
        mult = ONE/col_k[k];
        for(unsigned int i=k+1; i < m; i++) col_k[i] *= mult;
        
        /* row_i = row_i - [a(i,k)/a(k,k)] row_k, i=k+1, ..., m-1 */
        /* row k is the pivot row after swapping with row l.      */
        /* The computation is done one column at a time,          */
        /* column j=k+1, ..., n-1.                                */
        
        for (unsigned int j=k+1; j < n; j++) {
            
            col_j = a[j];
            a_kj = col_j[k];
            
            /* a(i,j) = a(i,j) - [a(i,k)/a(k,k)]*a(k,j)  */
            /* a_kj = a(k,j), col_k[i] = - a(i,k)/a(k,k) */
            
            if (a_kj != ZERO) {
                for (unsigned int i=k+1; i < m; i++)
                    col_j[i] -= a_kj * col_k[i];
            }
        }
    }
    
    /* return 0 to indicate success */
    
    return(0);
}

static void denseGETRS(double **a, unsigned int n, long int *p, double *b) {
    double *col_k, tmp;
    
    /* Permute b, based on pivot information in p */
    for (unsigned int k=0; k<n; k++) {
        long int pk = p[k];
        if(pk != k) {
            tmp = b[k];
            b[k] = b[pk];
            b[pk] = tmp;
        }
    }
    
    /* Solve Ly = b, store solution y in b */
    for (unsigned int k=0; k<n-1; k++) {
        col_k = a[k];
        for (unsigned int i=k+1; i<n; i++) b[i] -= col_k[i]*b[k];
    }
    
    /* Solve Ux = y, store solution x in b */
    for (unsigned int k = (unsigned int) n-1; k > 0; k--) {
        col_k = a[k];
        b[k] /= col_k[k];
        for (unsigned int i=0; i<k; i++) b[i] -= col_k[i]*b[k];
    }
    b[0] /= a[0][0];
}

static void denseCopy(double **a, double **b, long int m, long int n)
{
    long int i, j;
    double *a_col_j, *b_col_j;
    
    for (j=0; j < n; j++) {
        a_col_j = a[j];
        b_col_j = b[j];
        for (i=0; i < m; i++)
            b_col_j[i] = a_col_j[i];
    }
    
}

static void denseScale(double c, double **a, long int m, long int n)
{
    long int i, j;
    double *col_j;
    
    for (j=0; j < n; j++) {
        col_j = a[j];
        for (i=0; i < m; i++)
            col_j[i] *= c;
    }
}

static void denseAddIdentity(double **a, long int n)
{
    long int i;
    
    for (i=0; i < n; i++) a[i][i] += ONE;
}


static long int DenseGETRF(DlsMat A, long int *p)
{
    return(denseGETRF(A->cols, A->M, A->N, p));
}

static void DenseGETRS(DlsMat A, long int *p, double *b)
{
    denseGETRS(A->cols, A->N, p, b);
}

static void DenseCopy(DlsMat A, DlsMat B)
{
    denseCopy(A->cols, B->cols, A->M, A->N);
}

static void DenseScale(double c, DlsMat A)
{
    denseScale(c, A->cols, A->M, A->N);
}



/* Readability Replacements */

#define lmm       (cv_mem->cv_lmm)
#define f         (cv_mem->cv_f)
#define nst       (cv_mem->cv_nst)
#define tn        (cv_mem->cv_tn)
#define h         (cv_mem->cv_h)
#define gamma     (cv_mem->cv_gamma)
#define gammap    (cv_mem->cv_gammap)
#define gamrat    (cv_mem->cv_gamrat)
#define ewt       (cv_mem->cv_ewt)
#define linit     (cv_mem->cv_linit)
#define lsetup    (cv_mem->cv_lsetup)
//#define lsolve    (cv_mem->cv_lsolve)
#define lfree     (cv_mem->cv_lfree)
#define lmem      (cv_mem->cv_lmem)
#define vec_tmpl     (cv_mem->cv_tempv)
#define setupNonNull (cv_mem->cv_setupNonNull)

#define mtype     (cvdls_mem->d_type)
#define n         (cvdls_mem->d_n)
#define jacDQ     (cvdls_mem->d_jacDQ)
#define jac       (cvdls_mem->d_djac)
#define M         (cvdls_mem->d_M)
#define lpivots   (cvdls_mem->d_lpivots)
#define savedJ    (cvdls_mem->d_savedJ)
#define nstlj     (cvdls_mem->d_nstlj)
#define nje       (cvdls_mem->d_nje)
#define nfeDQ     (cvdls_mem->d_nfeDQ)
#define J_data    (cvdls_mem->d_J_data)
#define last_flag (cvdls_mem->d_last_flag)

                  
/*
 * -----------------------------------------------------------------
 * CVDense
 * -----------------------------------------------------------------
 * This routine initializes the memory record and sets various function
 * fields specific to the dense linear solver module.  CVDense first
 * calls the existing lfree routine if this is not NULL.  Then it sets
 * the cv_linit, cv_lsetup, cv_lsolve, cv_lfree fields in (*cvode_mem)
 * to be cvDenseInit, cvDenseSetup, cvDenseSolve, and cvDenseFree,
 * respectively.  It allocates memory for a structure of type
 * CVDlsMemRec and sets the cv_lmem field in (*cvode_mem) to the
 * address of this structure.  It sets setupNonNull in (*cvode_mem) to
 * TRUE, and the d_jac field to the default cvDlsDenseDQJac.
 * Finally, it allocates memory for M, savedJ, and lpivots.
 * The return value is SUCCESS = 0, or LMEM_FAIL = -1.
 *
 * NOTE: The dense linear solver assumes a serial implementation
 *       of the NVECTOR package. Therefore, CVDense will first 
 *       test for compatible a compatible N_Vector internal
 *       representation by checking that N_VGetArrayPointer and
 *       N_VSetArrayPointer exist.
 * -----------------------------------------------------------------
 */

int CVDense(void *cvode_mem, long int N)
{
  CVodeMem cv_mem;
  CVDlsMem cvdls_mem;

  /* Return immediately if cvode_mem is NULL */
  if (cvode_mem == NULL) {
    CVProcessError(NULL, CVDLS_MEM_NULL, "CVDENSE", "CVDense", MSGD_CVMEM_NULL);
    return(CVDLS_MEM_NULL);
  }
  cv_mem = (CVodeMem) cvode_mem;

  /* Test if the NVECTOR package is compatible with the DENSE solver */
  if (vec_tmpl->ops->nvgetarraypointer == NULL ||
      vec_tmpl->ops->nvsetarraypointer == NULL) {
    CVProcessError(cv_mem, CVDLS_ILL_INPUT, "CVDENSE", "CVDense", MSGD_BAD_NVECTOR);
    return(CVDLS_ILL_INPUT);
  }
    
  if (lfree !=NULL) lfree(cv_mem);

  /* Set four main function fields in cv_mem */
  linit  = cvDenseInit;
  lsetup = cvDenseSetup;
  cv_mem->cv_lsolve = cvDenseSolve;
  lfree  = cvDenseFree;

  /* Get memory for CVDlsMemRec */
  cvdls_mem = NULL;
  cvdls_mem = (CVDlsMem) malloc(sizeof(struct CVDlsMemRec));
  if (cvdls_mem == NULL) {
    CVProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    return(CVDLS_MEM_FAIL);
  }

  /* Set matrix type */
  mtype = SUNDIALS_DENSE;

  /* Initialize Jacobian-related data */
  jacDQ = TRUE;
  jac = NULL;
  J_data = NULL;

  last_flag = CVDLS_SUCCESS;

  setupNonNull = TRUE;

  /* Set problem dimension */
  n = N;

  /* Allocate memory for M, savedJ, and pivot array */

  M = NULL;
  M = NewDenseMat((unsigned int) N, (unsigned int) N);
  if (M == NULL) {
    CVProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  savedJ = NULL;
  savedJ = NewDenseMat((unsigned int) N, (unsigned int) N);
  if (savedJ == NULL) {
    CVProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }
  lpivots = NULL;
  lpivots = NewLintArray(N);
  if (lpivots == NULL) {
    CVProcessError(cv_mem, CVDLS_MEM_FAIL, "CVDENSE", "CVDense", MSGD_MEM_FAIL);
    DestroyMat(M);
    DestroyMat(savedJ);
    free(cvdls_mem); cvdls_mem = NULL;
    return(CVDLS_MEM_FAIL);
  }

  /* Attach linear solver memory to integrator memory */
  lmem = cvdls_mem;

  return(CVDLS_SUCCESS);
}

/*
 * -----------------------------------------------------------------
 * cvDenseInit
 * -----------------------------------------------------------------
 * This routine does remaining initializations specific to the dense
 * linear solver.
 * -----------------------------------------------------------------
 */

static int cvDenseInit(CVodeMem cv_mem)
{
  CVDlsMem cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  nje   = 0;
  nfeDQ = 0;
  nstlj = 0;

  /* Set Jacobian function and data, depending on jacDQ */
  if (jacDQ) {
    jac = cvDlsDenseDQJac;
    J_data = cv_mem;
  } else {
    J_data = cv_mem->cv_user_data;
  }

  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSetup
 * -----------------------------------------------------------------
 * This routine does the setup operations for the dense linear solver.
 * It makes a decision whether or not to call the Jacobian evaluation
 * routine based on various state variables, and if not it uses the 
 * saved copy.  In any case, it constructs the Newton matrix 
 * M = I - gamma*J, updates counters, and calls the dense LU 
 * factorization routine.
 * -----------------------------------------------------------------
 */

static int cvDenseSetup(CVodeMem cv_mem, int convfail, N_Vector ypred,
                        N_Vector fpred, int *jcurPtr, 
                        N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  int jbad, jok;
  double dgamma;
  long int ier;
  CVDlsMem cvdls_mem;
  int retval;

  cvdls_mem = (CVDlsMem) lmem;
 
  /* Use nst, gamma/gammap, and convfail to set J eval. flag jok */
 
  dgamma = fabs((gamma/gammap) - 1.0);
  jbad = (nst == 0) || (nst > nstlj + CVD_MSBJ) ||
         ((convfail == CV_FAIL_BAD_J) && (dgamma < CVD_DGMAX)) ||
         (convfail == CV_FAIL_OTHER);
  jok = !jbad;
 
  if (jok) {
    /* If jok = TRUE, use saved copy of J */
    *jcurPtr = FALSE;
    DenseCopy(savedJ, M);

  } else {

    /* If jok = FALSE, call jac routine for new J value */
    nje++;
    nstlj = nst;
    *jcurPtr = TRUE;
    SetToZero(M);

    retval = jac(n, tn, ypred, fpred, M, J_data, vtemp1, vtemp2, vtemp3);
    if (retval < 0) {
      CVProcessError(cv_mem, CVDLS_JACFUNC_UNRECVR, "CVDENSE", "cvDenseSetup", MSGD_JACFUNC_FAILED);
      last_flag = CVDLS_JACFUNC_UNRECVR;
      return(-1);
    }
    if (retval > 0) {
      last_flag = CVDLS_JACFUNC_RECVR;
      return(1);
    }

    DenseCopy(M, savedJ);

  }
  
  /* Scale and add I to get M = I - gamma*J */
  DenseScale(-gamma, M);
  AddIdentity(M);

  /* Do LU factorization of M */
  ier = DenseGETRF(M, lpivots); 

  /* Return 0 if the LU was complete; otherwise return 1 */
  last_flag = ier;
  if (ier > 0) return(1);
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseSolve
 * -----------------------------------------------------------------
 * This routine handles the solve operation for the dense linear solver
 * by calling the dense backsolve routine.  The returned value is 0.
 * -----------------------------------------------------------------
 */

static int cvDenseSolve(CVodeMem cv_mem, N_Vector b, N_Vector __attribute__((unused)) weight,
                        __attribute__((unused)) N_Vector ycur, __attribute__((unused)) N_Vector fcur) {
  CVDlsMem cvdls_mem;
  double *bd;

  cvdls_mem = (CVDlsMem) lmem;
  
  bd = N_VGetArrayPointer(b);

  DenseGETRS(M, lpivots, bd);

  /* If CV_BDF, scale the correction to account for change in gamma */
  if ((lmm == CV_BDF) && (gamrat != ONE)) {
    N_VScale(TWO/(ONE + gamrat), b, b);
  }
  
  last_flag = CVDLS_SUCCESS;
  return(0);
}

/*
 * -----------------------------------------------------------------
 * cvDenseFree
 * -----------------------------------------------------------------
 * This routine frees memory specific to the dense linear solver.
 * -----------------------------------------------------------------
 */

static void cvDenseFree(CVodeMem cv_mem)
{
  CVDlsMem  cvdls_mem;

  cvdls_mem = (CVDlsMem) lmem;
  
  DestroyMat(M);
  DestroyMat(savedJ);
  DestroyArray(lpivots);
  free(cvdls_mem);
  cv_mem->cv_lmem = NULL;
}


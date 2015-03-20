/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2010/12/01 22:10:38 $
 * ----------------------------------------------------------------- 
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * Common header file for the direct linear solvers in CVODE.
 * -----------------------------------------------------------------
 */

#ifndef _CVDLS_H
#define _CVDLS_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "sundials_direct.h"
#include "sundials_nvector.h"

/*
 * =================================================================
 *              C V D I R E C T     C O N S T A N T S
 * =================================================================
 */

/* 
 * -----------------------------------------------------------------
 * CVDLS return values 
 * -----------------------------------------------------------------
 */

#define CVDLS_SUCCESS           0
#define CVDLS_MEM_NULL         -1
#define CVDLS_LMEM_NULL        -2
#define CVDLS_ILL_INPUT        -3
#define CVDLS_MEM_FAIL         -4

/* Additional last_flag values */

#define CVDLS_JACFUNC_UNRECVR  -5
#define CVDLS_JACFUNC_RECVR    -6

#define TRUE 1
#define FALSE 0
    
/*
 * =================================================================
 *              F U N C T I O N   T Y P E S
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type: CVDlsDenseJacFn
 * -----------------------------------------------------------------
 *
 * A dense Jacobian approximation function Jac must be of type 
 * CVDlsDenseJacFn. Its parameters are:
 *
 * N   is the problem size.
 *
 * Jac is the dense matrix (of type DlsMat) that will be loaded
 *     by a CVDlsDenseJacFn with an approximation to the Jacobian 
 *     matrix J = (df_i/dy_j) at the point (t,y). 
 *
 * t   is the current value of the independent variable.
 *
 * y   is the current value of the dependent variable vector,
 *     namely the predicted value of y(t).
 *
 * fy  is the vector f(t,y).
 *
 * user_data is a pointer to user data - the same as the user_data
 *     parameter passed to CVodeSetFdata.
 *
 * tmp1, tmp2, and tmp3 are pointers to memory allocated for
 * vectors of length N which can be used by a CVDlsDenseJacFn
 * as temporary storage or work space.
 *
 * A CVDlsDenseJacFn should return 0 if successful, a positive 
 * value if a recoverable error occurred, and a negative value if 
 * an unrecoverable error occurred.
 *
 * -----------------------------------------------------------------
 *
 * NOTE: The following are two efficient ways to load a dense Jac:         
 * (1) (with macros - no explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = DENSE_COL(Jac,j);                                 
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * (2) (without macros - explicit data structure references)      
 *     for (j=0; j < Neq; j++) {                                  
 *       col_j = (Jac->data)[j];                                   
 *       for (i=0; i < Neq; i++) {                                
 *         generate J_ij = the (i,j)th Jacobian element           
 *         col_j[i] = J_ij;                                       
 *       }                                                        
 *     }                                                          
 * A third way, using the DENSE_ELEM(A,i,j) macro, is much less   
 * efficient in general.  It is only appropriate for use in small 
 * problems in which efficiency of access is NOT a major concern. 
 *                                                                
 * NOTE: If the user's Jacobian routine needs other quantities,   
 *     they are accessible as follows: hcur (the current stepsize)
 *     and ewt (the error weight vector) are accessible through   
 *     CVodeGetCurrentStep and CVodeGetErrWeights, respectively 
 *     (see cvode.h). The unit roundoff is available as 
 *     UNIT_ROUNDOFF defined in sundials_types.h.
 *
 * -----------------------------------------------------------------
 */
  
  
typedef int (*CVDlsDenseJacFn)(long int N, double t,
			       N_Vector y, N_Vector fy, 
			       DlsMat Jac, void *user_data,
			       N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  


/*
 * =================================================================
 *            E X P O R T E D    F U N C T I O N S 
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Optional inputs to the CVDLS linear solver
 * -----------------------------------------------------------------
 *
 * CVDlsSetDenseJacFn specifies the dense Jacobian approximation
 * routine to be used for a direct dense linear solver.
 *
 * By default, a difference quotient approximation, supplied with
 * the solver is used.
 *
 * The return value is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODE memory was NULL
 *    CVDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDlsSetDenseJacFn(void *cvode_mem, CVDlsDenseJacFn jac);

/*
 * -----------------------------------------------------------------
 * Optional outputs from the CVDLS linear solver
 * -----------------------------------------------------------------
 *
 * CVDlsGetNumJacEvals returns the number of calls made to the
 *                     Jacobian evaluation routine jac.
 * CVDlsGetNumRhsEvals returns the number of calls to the user
 *                     f routine due to finite difference Jacobian
 *                     evaluation.
 * CVDlsGetLastFlag    returns the last error flag set by any of
 *                     the CVDLS interface functions.
 *
 * The return value of CVDlsGet* is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the CVODE memory was NULL
 *    CVDLS_LMEM_NULL if the linear solver memory was NULL
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDlsGetNumJacEvals(void *cvode_mem, long int *njevals);
SUNDIALS_EXPORT int CVDlsGetNumRhsEvals(void *cvode_mem, long int *nfevalsLS);
SUNDIALS_EXPORT int CVDlsGetLastFlag(void *cvode_mem, long int *flag);

/*
 * -----------------------------------------------------------------
 * The following function returns the name of the constant 
 * associated with a CVDLS return flag
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT char *CVDlsGetReturnFlagName(long int flag);


#ifdef __cplusplus
}
#endif

#endif

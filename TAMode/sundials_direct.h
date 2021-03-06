/*
 * -----------------------------------------------------------------
 * $Revision: 1.6 $
 * $Date: 2010/12/22 22:18:49 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This header file contains definitions and declarations for use by
 * generic direct linear solvers for Ax = b. It defines types for
 * dense and banded matrices and corresponding accessor macros.
 * -----------------------------------------------------------------
 */

#ifndef _SUNDIALS_DIRECT_H
#define _SUNDIALS_DIRECT_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif
    
    
#include <float.h>

/*
 * ==================================================================
 * Type definitions
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * Type : DlsMat
 * -----------------------------------------------------------------
 * The type DlsMat is defined to be a pointer to a structure
 * with various sizes, a data field, and an array of pointers to
 * the columns which defines a dense or band matrix for use in 
 * direct linear solvers. The M and N fields indicates the number 
 * of rows and columns, respectively. The data field is a one 
 * dimensional array used for component storage. The cols field 
 * stores the pointers in data for the beginning of each column.
 * -----------------------------------------------------------------
 * For DENSE matrices, the relevant fields in DlsMat are:
 *    M     - number of rows
 *    N     - number of columns
 *    ldim  - leading dimension (ldim >= M)
 *    data  - pointer to a contiguous block of double variables
 *    ldata - length of the data array =ldim*N
 *    cols  - array of pointers. cols[j] points to the first element 
 *            of the j-th column of the matrix in the array data.
 *
 * The elements of a dense matrix are stored columnwise (i.e columns 
 * are stored one on top of the other in memory). 
 * If A is of type DlsMat, then the (i,j)th element of A (with 
 * 0 <= i < M and 0 <= j < N) is given by (A->data)[j*n+i]. 
 *
 * The DENSE_COL and DENSE_ELEM macros below allow a user to access 
 * efficiently individual matrix elements without writing out explicit 
 * data structure references and without knowing too much about the 
 * underlying element storage. The only storage assumption needed is 
 * that elements are stored columnwise and that a pointer to the 
 * jth column of elements can be obtained via the DENSE_COL macro.
 * -----------------------------------------------------------------
 * For BAND matrices, the relevant fields in DlsMat are:
 *    type  = SUNDIALS_BAND
 *    M     - number of rows
 *    N     - number of columns
 *    mu    - upper bandwidth, 0 <= mu <= min(M,N)
 *    ml    - lower bandwidth, 0 <= ml <= min(M,N)
 *    s_mu  - storage upper bandwidth, mu <= s_mu <= N-1.
 *            The dgbtrf routine writes the LU factors into the storage 
 *            for A. The upper triangular factor U, however, may have 
 *            an upper bandwidth as big as MIN(N-1,mu+ml) because of 
 *            partial pivoting. The s_mu field holds the upper 
 *            bandwidth allocated for A.
 *    ldim  - leading dimension (ldim >= s_mu)
 *    data  - pointer to a contiguous block of double variables
 *    ldata - length of the data array =ldim*(s_mu+ml+1)
 *    cols  - array of pointers. cols[j] points to the first element 
 *            of the j-th column of the matrix in the array data.
 *
 * The BAND_COL, BAND_COL_ELEM, and BAND_ELEM macros below allow a 
 * user to access individual matrix elements without writing out 
 * explicit data structure references and without knowing too much 
 * about the underlying element storage. The only storage assumption 
 * needed is that elements are stored columnwise and that a pointer 
 * into the jth column of elements can be obtained via the BAND_COL 
 * macro. The BAND_COL_ELEM macro selects an element from a column
 * which has already been isolated via BAND_COL. The macro 
 * BAND_COL_ELEM allows the user to avoid the translation 
 * from the matrix location (i,j) to the index in the array returned 
 * by BAND_COL at which the (i,j)th element is stored. 
 * -----------------------------------------------------------------
 */

typedef struct _DlsMat {
  int type;
  unsigned int M;
  unsigned int N;
  long int ldim;
  long int mu;
  long int ml;
  long int s_mu;
  double *data;
  long int ldata;
  double **cols;
} *DlsMat;

/*
 * ==================================================================
 * Data accessor macros
 * ==================================================================
 */

/*
 * -----------------------------------------------------------------
 * DENSE_COL and DENSE_ELEM
 * -----------------------------------------------------------------
 *
 * DENSE_COL(A,j) references the jth column of the M-by-N dense
 * matrix A, 0 <= j < N. The type of the expression DENSE_COL(A,j) 
 * is (double *). After the assignment in the usage above, col_j 
 * may be treated as an array indexed from 0 to M-1. The (i,j)-th 
 * element of A is thus referenced by col_j[i].
 *
 * DENSE_ELEM(A,i,j) references the (i,j)th element of the dense 
 * M-by-N matrix A, 0 <= i < M ; 0 <= j < N.
 *
 * -----------------------------------------------------------------
 */

#define DENSE_COL(A,j) ((A->cols)[j])
#define DENSE_ELEM(A,i,j) ((A->cols)[j][i])

/*
 * -----------------------------------------------------------------
 * BAND_COL, BAND_COL_ELEM, and BAND_ELEM
 * -----------------------------------------------------------------
 *  
 * BAND_COL(A,j) references the diagonal element of the jth column 
 * of the N by N band matrix A, 0 <= j <= N-1. The type of the 
 * expression BAND_COL(A,j) is double *. The pointer returned by 
 * the call BAND_COL(A,j) can be treated as an array which is 
 * indexed from -(A->mu) to (A->ml).
 * 
 * BAND_COL_ELEM references the (i,j)th entry of the band matrix A 
 * when used in conjunction with BAND_COL. The index (i,j) should 
 * satisfy j-(A->mu) <= i <= j+(A->ml).
 *
 * BAND_ELEM(A,i,j) references the (i,j)th element of the M-by-N 
 * band matrix A, where 0 <= i,j <= N-1. The location (i,j) should 
 * further satisfy j-(A->mu) <= i <= j+(A->ml). 
 *
 * -----------------------------------------------------------------
 */
 
#define BAND_COL(A,j) (((A->cols)[j])+(A->s_mu))
#define BAND_COL_ELEM(col_j,i,j) (col_j[(i)-(j)])
#define BAND_ELEM(A,i,j) ((A->cols)[j][(i)-(j)+(A->s_mu)])

DlsMat NewDenseMat(unsigned int M, unsigned int N);
void DestroyMat(DlsMat A);
int *NewIntArray(int N);
long int *NewLintArray(unsigned long int N);
double *NewRealArray(long int N);
void DestroyArray(void *p);
void AddIdentity(DlsMat A);
void SetToZero(DlsMat A);
void PrintMat(DlsMat A);

double **newDenseMat(unsigned int m, unsigned int n);
void destroyMat(double **a);
int *newIntArray(int n);
long int *newLintArray(long int n);
double *newRealArray(long int m);
void destroyArray(void *v);


#ifdef __cplusplus
}
#endif

#endif

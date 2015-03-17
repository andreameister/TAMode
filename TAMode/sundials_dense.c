/*
 * -----------------------------------------------------------------
 * $Revision: 1.7 $
 * $Date: 2010/12/01 22:46:56 $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic package of dense
 * matrix operations.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>

#include "sundials_dense.h"
#include "sundials_nvector.h"

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)

/*
 * -----------------------------------------------------
 * Functions working on DlsMat
 * -----------------------------------------------------
 */

long int DenseGETRF(DlsMat A, long int *p)
{
    return(denseGETRF(A->cols, A->M, A->N, p));
}

void DenseGETRS(DlsMat A, long int *p, double *b)
{
    denseGETRS(A->cols, A->N, p, b);
}

long int DensePOTRF(DlsMat A)
{
    return(densePOTRF(A->cols, A->M));
}

void DensePOTRS(DlsMat A, double *b)
{
    densePOTRS(A->cols, A->M, b);
}

int DenseGEQRF(DlsMat A, double *beta, double *wrk)
{
    return(denseGEQRF(A->cols, A->M, A->N, beta, wrk));
}

void DenseCopy(DlsMat A, DlsMat B)
{
    denseCopy(A->cols, B->cols, A->M, A->N);
}

void DenseScale(double c, DlsMat A)
{
    denseScale(c, A->cols, A->M, A->N);
}

long int denseGETRF(double **a, unsigned int m, unsigned int n, long int *p) {
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

void denseGETRS(double **a, unsigned int n, long int *p, double *b) {
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

/*
 * Cholesky decomposition of a symmetric positive-definite matrix
 * A = C^T*C: gaxpy version.
 * Only the lower triangle of A is accessed and it is overwritten with
 * the lower triangle of C.
 */
long int densePOTRF(double **a, long int m)
{
    double *a_col_j, *a_col_k;
    double a_diag;
    long int i, j, k;
    
    for (j=0; j<m; j++) {
        
        a_col_j = a[j];
        
        if (j>0) {
            for(i=j; i<m; i++) {
                for(k=0;k<j;k++) {
                    a_col_k = a[k];
                    a_col_j[i] -= a_col_k[i]*a_col_k[j];
                }
            }
        }
        
        a_diag = a_col_j[j];
        if (a_diag <= ZERO) return(j+1);
        a_diag = RSqrt(a_diag);
        
        for(i=j; i<m; i++) a_col_j[i] /= a_diag;
        
    }
    
    return(0);
}

/*
 * Solution of Ax=b, with A s.p.d., based on the Cholesky decomposition
 * obtained with denPOTRF.; A = C*C^T, C lower triangular
 *
 */
void densePOTRS(double **a, long int m, double *b)
{
    double *col_j, *col_i;
    long int i, j;
    
    /* Solve C y = b, forward substitution - column version.
     Store solution y in b */
    for (j=0; j < m-1; j++) {
        col_j = a[j];
        b[j] /= col_j[j];
        for (i=j+1; i < m; i++)
            b[i] -= b[j]*col_j[i];
    }
    col_j = a[m-1];
    b[m-1] /= col_j[m-1];
    
    /* Solve C^T x = y, backward substitution - row version.
     Store solution x in b */
    col_j = a[m-1];
    b[m-1] /= col_j[m-1];
    for (i=m-2; i>=0; i--) {
        col_i = a[i];
        for (j=i+1; j<m; j++)
            b[i] -= col_i[j]*b[j];
        b[i] /= col_i[i];
    }
    
}

/*
 * QR factorization of a rectangular matrix A of size m by n (m >= n)
 * using Householder reflections.
 *
 * On exit, the elements on and above the diagonal of A contain the n by n
 * upper triangular matrix R; the elements below the diagonal, with the array beta,
 * represent the orthogonal matrix Q as a product of elementary reflectors .
 *
 * v (of length m) must be provided as workspace.
 *
 */

int denseGEQRF(double **a, long int m, long int n, double *beta, double *v)
{
    double ajj, s, mu, v1, v1_2;
    double *col_j, *col_k;
    long int i, j, k;
    
    /* For each column...*/
    for(j=0; j<n; j++) {
        
        col_j = a[j];
        
        ajj = col_j[j];
        
        /* Compute the j-th Householder vector (of length m-j) */
        v[0] = ONE;
        s = ZERO;
        for(i=1; i<m-j; i++) {
            v[i] = col_j[i+j];
            s += v[i]*v[i];
        }
        
        if(s != ZERO) {
            mu = RSqrt(ajj*ajj+s);
            v1 = (ajj <= ZERO) ? ajj-mu : -s/(ajj+mu);
            v1_2 = v1*v1;
            beta[j] = TWO * v1_2 / (s + v1_2);
            for(i=1; i<m-j; i++) v[i] /= v1;
        } else {
            beta[j] = ZERO;
        }
        
        /* Update upper triangle of A (load R) */
        for(k=j; k<n; k++) {
            col_k = a[k];
            s = ZERO;
            for(i=0; i<m-j; i++) s += col_k[i+j]*v[i];
            s *= beta[j];
            for(i=0; i<m-j; i++) col_k[i+j] -= s*v[i];
        }
        
        /* Update A (load Householder vector) */
        if(j<m-1) {
            for(i=1; i<m-j; i++) col_j[i+j] = v[i];
        }
        
    }
    
    
    return(0);
}

void denseCopy(double **a, double **b, long int m, long int n)
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

void denseScale(double c, double **a, long int m, long int n)
{
    long int i, j;
    double *col_j;
    
    for (j=0; j < n; j++) {
        col_j = a[j];
        for (i=0; i < m; i++)
            col_j[i] *= c;
    }
}

void denseAddIdentity(double **a, long int n)
{
    long int i;
    
    for (i=0; i < n; i++) a[i][i] += ONE;
}

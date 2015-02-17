/* cobyla : contrained optimization by linear approximation */

/*
 * Copyright (c) 1992, Michael J. D. Powell (M.J.D.Powell@damtp.cam.ac.uk)
 * Copyright (c) 2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C version of COBYLA2, a contrained optimization by linear
 * approximation package developed by Michael J. D. Powell in Fortran.
 * 
 * The original source code can be found at :
 * http://plato.la.asu.edu/topics/problems/nlores.html
 */

/* $Jeannot: cobyla.h,v 1.10 2004/04/18 09:51:37 js Exp $ */

#ifndef _COBYLA_
#define _COBYLA_

/* A preconditioner, which preconditions v at x to return vpre.
 (The meaning of "preconditioning" is algorithm-dependent.) */
typedef void (*nlopt_precond)(unsigned n, const double *x, const double *v,
double *vpre, void *data);

typedef double (*nlopt_func)(unsigned n, const double *x,
double *gradient, /* NULL if not needed */
void *func_data);

typedef void (*nlopt_mfunc)(unsigned m, double *result,
unsigned n, const double *x,
double *gradient, /* NULL if not needed */
void *func_data);



/* stopping criteria */
typedef struct {
    unsigned n;
    double minf_max;
    double ftol_rel;
    double ftol_abs;
    double xtol_rel;
    const double *xtol_abs;
    int nevals, maxeval;
    double maxtime, start;
    int *force_stop;
} nlopt_stopping;

typedef struct {
    unsigned m; /* dimensional of constraint: mf maps R^n -> R^m */
    nlopt_func f; /* one-dimensional constraint, requires m == 1 */
    nlopt_mfunc mf;
    nlopt_precond pre; /* preconditioner for f (NULL if none or if mf) */
    void *f_data;
    double *tol;
} nlopt_constraint;

typedef enum {
    NLOPT_FAILURE = -1, /* generic failure code */
    NLOPT_INVALID_ARGS = -2,
    NLOPT_OUT_OF_MEMORY = -3,
    NLOPT_ROUNDOFF_LIMITED = -4,
    NLOPT_FORCED_STOP = -5,
    NLOPT_SUCCESS = 1, /* generic success code */
    NLOPT_STOPVAL_REACHED = 2,
    NLOPT_FTOL_REACHED = 3,
    NLOPT_XTOL_REACHED = 4,
    NLOPT_MAXEVAL_REACHED = 5,
    NLOPT_MAXTIME_REACHED = 6
} nlopt_result;

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

/* NLopt-style interface function */
nlopt_result cobyla_minimize(unsigned n, nlopt_func f, void *f_data,
                             unsigned m, nlopt_constraint *fc,
                             unsigned p, nlopt_constraint *h,
                             const double *lb, const double *ub, /* bounds */
                             double *x, /* in: initial guess, out: minimizer */
                             double *minf,
                             nlopt_stopping *stop,
                             const double *dx);

#ifdef __cplusplus
}
#endif

#endif /* _COBYLA_ */

/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2007/04/06 20:33:30 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Radu Serban and Aaron Collier @ LLNL                               
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic NVECTOR package.
 * It contains the implementation of the N_Vector operations listed
 * in nvector.h.
 * -----------------------------------------------------------------
 */

#include <stdlib.h>

#include "sundials_nvector.h"

#define ZERO (0.0)
#define ONE  (1.0)

/*
 * -----------------------------------------------------------------
 * Functions in the 'ops' structure
 * -----------------------------------------------------------------
 */

double RPowerI(double base, int exponent)
{
    int i, expt;
    double prod;
    
    prod = ONE;
    expt = abs(exponent);
    for(i = 1; i <= expt; i++) prod *= base;
    if (exponent < 0) prod = ONE/prod;
    return(prod);
}

double RPowerR(double base, double exponent) {
    if (base <= ZERO) return(ZERO);
    
    return(pow(base, exponent));
}

double RSqrt(double x) {
    if (x <= ZERO) return(ZERO);
    return(sqrt(x));
}

double RExp(double x)
{
    return(exp(x));
}


N_Vector N_VClone(N_Vector w)
{
  N_Vector v = NULL;
  v = w->ops->nvclone(w);
  return(v);
}

N_Vector N_VCloneEmpty(N_Vector w)
{
  N_Vector v = NULL;
  v = w->ops->nvcloneempty(w);
  return(v);
}

void N_VDestroy(N_Vector v)
{
  if (v==NULL) return;
  v->ops->nvdestroy(v);
  return;
}

void N_VSpace(N_Vector v, long int *lrw, long int *liw)
{
  v->ops->nvspace(v, lrw, liw);
  return;
}

double *N_VGetArrayPointer(N_Vector v)
{
  return((double *) v->ops->nvgetarraypointer(v));
}

void N_VSetArrayPointer(double *v_data, N_Vector v)
{
  v->ops->nvsetarraypointer(v_data, v);
  return;
}

void N_VLinearSum(double a, N_Vector x, double b, N_Vector y, N_Vector z)
{
  z->ops->nvlinearsum(a, x, b, y, z);
  return;
}

void N_VConst(double c, N_Vector z)
{
  z->ops->nvconst(c, z);
  return;
}

void N_VProd(N_Vector x, N_Vector y, N_Vector z)
{
  z->ops->nvprod(x, y, z);
  return;
}

void N_VDiv(N_Vector x, N_Vector y, N_Vector z)
{
  z->ops->nvdiv(x, y, z);
  return;
}

void N_VScale(double c, N_Vector x, N_Vector z) 
{
  z->ops->nvscale(c, x, z);
  return;
}

void N_VAbs(N_Vector x, N_Vector z)
{
  z->ops->nvabs(x, z);
  return;
}

void N_VInv(N_Vector x, N_Vector z)
{
  z->ops->nvinv(x, z);
  return;
}

void N_VAddConst(N_Vector x, double b, N_Vector z)
{
  z->ops->nvaddconst(x, b, z);
  return;
}

double N_VDotProd(N_Vector x, N_Vector y)
{
  return((double) y->ops->nvdotprod(x, y));
}

double N_VMaxNorm(N_Vector x)
{
  return((double) x->ops->nvmaxnorm(x));
}

double N_VWrmsNorm(N_Vector x, N_Vector w)
{
  return((double) x->ops->nvwrmsnorm(x, w));
}

double N_VMin(N_Vector x)
{
  return((double) x->ops->nvmin(x));
}

void N_VCompare(double c, N_Vector x, N_Vector z)
{
  z->ops->nvcompare(c, x, z);
  return;
}

/*
 * -----------------------------------------------------------------
 * Additional functions exported by the generic NVECTOR:
 *   N_VCloneEmptyVectorArray
 *   N_VCloneVectorArray
 *   N_VDestroyVectorArray
 * -----------------------------------------------------------------
 */

N_Vector *N_VCloneEmptyVectorArray(int count, N_Vector w)
{
  N_Vector *vs = NULL;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc((size_t) count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VCloneEmpty(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

N_Vector *N_VCloneVectorArray(int count, N_Vector w)
{
  N_Vector *vs = NULL;
  int j;

  if (count <= 0) return(NULL);

  vs = (N_Vector *) malloc((size_t) count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = N_VClone(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

void N_VDestroyVectorArray(N_Vector *vs, int count)
{
  int j;

  if (vs==NULL) return;

  for (j = 0; j < count; j++) N_VDestroy(vs[j]);

  free(vs); vs = NULL;

  return;
}

/*
 * -----------------------------------------------------------------
 * $Revision: 1.1 $
 * $Date: 2006/07/05 15:32:37 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh, Radu Serban,
 *                and Aaron Collier @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a serial implementation
 * of the NVECTOR package.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "nvector_serial.h"

#define ZERO   (0.0)
#define HALF   (0.5)
#define ONE    (1.0)
#define ONEPT5 (1.5)
#define FALSE 0
#define TRUE 1

/* Private function prototypes */
/* z=x */
static void VCopy_Serial(N_Vector x, N_Vector z);
/* z=x+y */
static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z);
/* z=x-y */
static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z);
/* z=-x */
static void VNeg_Serial(N_Vector x, N_Vector z);
/* z=c(x-y) */
static void VScaleDiff_Serial(double c, N_Vector x, N_Vector y, N_Vector z); 
/* z=ax+y */
static void VLin1_Serial(double a, N_Vector x, N_Vector y, N_Vector z);
/* z=ax-y */
static void VLin2_Serial(double a, N_Vector x, N_Vector y, N_Vector z);
/* y <- ax+y */
static void Vaxpy_Serial(double a, N_Vector x, N_Vector y);
/* x <- ax */
static void VScaleBy_Serial(double a, N_Vector x);

/*
 * -----------------------------------------------------------------
 * exported functions
 * -----------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------
 * Function to create a new empty serial vector 
 */

N_Vector N_VNewEmpty_Serial(long int length)
{
  N_Vector v;
  N_Vector_Ops ops;
  N_VectorContent_Serial content;

  /* Create vector */
  v = NULL;
  v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);
  
  /* Create vector operation structure */
  ops = NULL;
  ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }

  ops->nvclone           = N_VClone_Serial;
  ops->nvcloneempty      = N_VCloneEmpty_Serial;
  ops->nvdestroy         = N_VDestroy_Serial;
  ops->nvspace           = N_VSpace_Serial;
  ops->nvgetarraypointer = N_VGetArrayPointer_Serial;
  ops->nvsetarraypointer = N_VSetArrayPointer_Serial;
  ops->nvlinearsum       = N_VLinearSum_Serial;
  ops->nvconst           = N_VConst_Serial;
  ops->nvprod            = N_VProd_Serial;
  ops->nvdiv             = N_VDiv_Serial;
  ops->nvscale           = N_VScale_Serial;
  ops->nvabs             = N_VAbs_Serial;
  ops->nvinv             = N_VInv_Serial;
  ops->nvaddconst        = N_VAddConst_Serial;
  ops->nvdotprod         = N_VDotProd_Serial;
  ops->nvmaxnorm         = N_VMaxNorm_Serial;
  ops->nvwrmsnorm        = N_VWrmsNorm_Serial;
  ops->nvmin             = N_VMin_Serial;
  ops->nvcompare         = N_VCompare_Serial;

  /* Create content */
  content = NULL;
  content = (N_VectorContent_Serial) malloc(sizeof(struct _N_VectorContent_Serial));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length   = length;
  content->own_data = FALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a new serial vector 
 */

N_Vector N_VNew_Serial(long int length)
{
  N_Vector v;
  double *data;

  v = NULL;
  v = N_VNewEmpty_Serial(length);
  if (v == NULL) return(NULL);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (double *) malloc((size_t) length * sizeof(double));
    if(data == NULL) { N_VDestroy_Serial(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_S(v) = TRUE;
    NV_DATA_S(v)     = data;

  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create a serial N_Vector with user data component 
 */

N_Vector N_VMake_Serial(long int length, double *v_data)
{
  N_Vector v;

  v = NULL;
  v = N_VNewEmpty_Serial(length);
  if (v == NULL) return(NULL);

  if (length > 0) {
    /* Attach data */
    NV_OWN_DATA_S(v) = FALSE;
    NV_DATA_S(v)     = v_data;
  }

  return(v);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new serial vectors. 
 */

N_Vector *N_VCloneVectorArray_Serial(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc((size_t) count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VClone_Serial(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Serial(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to create an array of new serial vectors with NULL data array. 
 */

N_Vector *N_VCloneVectorArrayEmpty_Serial(int count, N_Vector w)
{
  N_Vector *vs;
  int j;

  if (count <= 0) return(NULL);

  vs = NULL;
  vs = (N_Vector *) malloc((size_t) count * sizeof(N_Vector));
  if(vs == NULL) return(NULL);

  for (j = 0; j < count; j++) {
    vs[j] = NULL;
    vs[j] = N_VCloneEmpty_Serial(w);
    if (vs[j] == NULL) {
      N_VDestroyVectorArray_Serial(vs, j-1);
      return(NULL);
    }
  }

  return(vs);
}

/* ----------------------------------------------------------------------------
 * Function to free an array created with N_VCloneVectorArray_Serial
 */

void N_VDestroyVectorArray_Serial(N_Vector *vs, int count)
{
  int j;

  for (j = 0; j < count; j++) N_VDestroy_Serial(vs[j]);

  free(vs); vs = NULL;

  return;
}

/* ----------------------------------------------------------------------------
 * Function to print the a serial vector 
 */
 
void N_VPrint_Serial(N_Vector x)
{
  long int i, N;
  double *xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++) {
    printf("%11.8lg\n", xd[i]);
  }
  printf("\n");

  return;
}

/*
 * -----------------------------------------------------------------
 * implementation of vector operations
 * -----------------------------------------------------------------
 */

N_Vector N_VCloneEmpty_Serial(N_Vector w) {
  if (w == NULL) return(NULL);

  /* Create vector */
  N_Vector v = (N_Vector) malloc(sizeof *v);
  if (v == NULL) return(NULL);

  /* Create vector operation structure */
  N_Vector_Ops ops = (N_Vector_Ops) malloc(sizeof(struct _generic_N_Vector_Ops));
  if (ops == NULL) { free(v); return(NULL); }
  
  ops->nvclone           = w->ops->nvclone;
  ops->nvcloneempty      = w->ops->nvcloneempty;
  ops->nvdestroy         = w->ops->nvdestroy;
  ops->nvspace           = w->ops->nvspace;
  ops->nvgetarraypointer = w->ops->nvgetarraypointer;
  ops->nvsetarraypointer = w->ops->nvsetarraypointer;
  ops->nvlinearsum       = w->ops->nvlinearsum;
  ops->nvconst           = w->ops->nvconst;  
  ops->nvprod            = w->ops->nvprod;   
  ops->nvdiv             = w->ops->nvdiv;
  ops->nvscale           = w->ops->nvscale; 
  ops->nvabs             = w->ops->nvabs;
  ops->nvinv             = w->ops->nvinv;
  ops->nvaddconst        = w->ops->nvaddconst;
  ops->nvdotprod         = w->ops->nvdotprod;
  ops->nvmaxnorm         = w->ops->nvmaxnorm;
  ops->nvwrmsnorm        = w->ops->nvwrmsnorm;
  ops->nvmin             = w->ops->nvmin;
  ops->nvcompare         = w->ops->nvcompare;

  /* Create content */
  N_VectorContent_Serial content = (N_VectorContent_Serial) malloc(sizeof(struct _N_VectorContent_Serial));
  if (content == NULL) { free(ops); free(v); return(NULL); }

  content->length   = NV_LENGTH_S(w);
  content->own_data = FALSE;
  content->data     = NULL;

  /* Attach content and ops */
  v->content = content;
  v->ops     = ops;

  return(v);
}

N_Vector N_VClone_Serial(N_Vector w)
{
  N_Vector v;
  double *data;
  long int length;

  v = NULL;
  v = N_VCloneEmpty_Serial(w);
  if (v == NULL) return(NULL);

  length = NV_LENGTH_S(w);

  /* Create data */
  if (length > 0) {

    /* Allocate memory */
    data = NULL;
    data = (double *) malloc((size_t) length * sizeof(double));
    if(data == NULL) { N_VDestroy_Serial(v); return(NULL); }

    /* Attach data */
    NV_OWN_DATA_S(v) = TRUE;
    NV_DATA_S(v)     = data;

  }

  return(v);
}

void N_VDestroy_Serial(N_Vector v)
{
  if (NV_OWN_DATA_S(v) == TRUE) {
    free(NV_DATA_S(v));
    NV_DATA_S(v) = NULL;
  }
  free(v->content); v->content = NULL;
  free(v->ops); v->ops = NULL;
  free(v); v = NULL;

  return;
}

void N_VSpace_Serial(N_Vector v, long int *lrw, long int *liw)
{
  *lrw = NV_LENGTH_S(v);
  *liw = 1;

  return;
}

double *N_VGetArrayPointer_Serial(N_Vector v)
{
  return((double *) NV_DATA_S(v));
}

void N_VSetArrayPointer_Serial(double *v_data, N_Vector v)
{
  if (NV_LENGTH_S(v) > 0) NV_DATA_S(v) = v_data;

  return;
}

void N_VLinearSum_Serial(double a, N_Vector x, double b, N_Vector y, N_Vector z)
{
    if ((b == 1.0) && (z == y)) {    /* BLAS usage: axpy y <- ax+y */
        Vaxpy_Serial(a,x,y);
        return;
    }
    
    if ((a == 1.0) && (z == x)) {    /* BLAS usage: axpy x <- by+x */
        Vaxpy_Serial(b,y,x);
        return;
    }
    
    /* Case: a == b == 1.0 */
    
    if ((a == 1.0) && (b == 1.0)) {
        VSum_Serial(x, y, z);
        return;
    }
    
    N_Vector v1, v2;
    int test;
    /* Cases: (1) a == 1.0, b = -1.0, (2) a == -1.0, b == 1.0 */
    
    if ((test = ((a == ONE) && (b == -ONE))) || ((a == -ONE) && (b == ONE))) {
        v1 = test ? y : x;
        v2 = test ? x : y;
        VDiff_Serial(v2, v1, z);
        return;
    }
    
    /* Cases: (1) a == 1.0, b == other or 0.0, (2) a == other or 0.0, b == 1.0 */
    /* if a or b is 0.0, then user should have called N_VScale */
    
    if ((test = (a == ONE)) || (b == ONE)) {
        v1 = test ? y : x;
        v2 = test ? x : y;
        VLin1_Serial(test ? b : a, v1, v2, z);
        return;
    }
    
    /* Cases: (1) a == -1.0, b != 1.0, (2) a != 1.0, b == -1.0 */
    
    if ((test = (a == -ONE)) || (b == -ONE)) {
        v1 = test ? y : x;
        v2 = test ? x : y;
        VLin2_Serial(test ? b : a, v1, v2, z);
        return;
    }
    
    /* Case: a == b */
    /* catches case both a and b are 0.0 - user should have called N_VConst */
    
    if (a == b) {
        const long int N  = NV_LENGTH_S(x);
        const double *xd = NV_DATA_S(x);
        const double *yd = NV_DATA_S(y);
        double *zd = NV_DATA_S(z);
        
        for (size_t i = 0; i < N; i++)
            zd[i] = a*(xd[i]+yd[i]);
        
        return;
    }
    
    /* Case: a == -b */
    
    if (a == -b) {
        VScaleDiff_Serial(a, x, y, z);
        return;
    }
    
    /* Do all cases not handled above:
     (1) a == other, b == 0.0 - user should have called N_VScale
     (2) a == 0.0, b == other - user should have called N_VScale
     (3) a,b == other, a !=b, a != -b */
    
    const long int N  = NV_LENGTH_S(x);
    const double *xd = NV_DATA_S(x);
    const double *yd = NV_DATA_S(y);
    double *zd = NV_DATA_S(z);
    
    for (size_t i = 0; i < N; i++)
        zd[i] = (a*xd[i])+(b*yd[i]);
    
    return;
}

void N_VConst_Serial(double c, N_Vector z)
{
  long int i, N;
  double *zd;

  zd = NULL;

  N  = NV_LENGTH_S(z);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) zd[i] = c;

  return;
}

void N_VProd_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  double *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]*yd[i];

  return;
}

void N_VDiv_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  double *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]/yd[i];

  return;
}

void N_VScale_Serial(const double c, N_Vector x, N_Vector z)
{
  if (z == x) {  /* BLAS usage: scale x <- cx */
    VScaleBy_Serial(c, x);
    return;
  }

  if (c == ONE) {
    VCopy_Serial(x, z);
  } else if (c == -ONE) {
    VNeg_Serial(x, z);
  } else {
    const long int N  = NV_LENGTH_S(x);
    const double * const xd = NV_DATA_S(x);
    double * const zd = NV_DATA_S(z);
    for (size_t i = 0; i < N; i++)
      zd[i] = c*xd[i];
  }

  return;
}

void N_VAbs_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  double *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++)
    zd[i] = fabs(xd[i]);

  return;
}

void N_VInv_Serial(N_Vector x, N_Vector z)
{
  double *xd, *zd;

  long N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (size_t i = 0; i < N; i++)
    zd[i] = ONE/xd[i];

  return;
}

void N_VAddConst_Serial(N_Vector x, double b, N_Vector z)
{
  long int i, N;
  double *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) 
    zd[i] = xd[i]+b;

  return;
}

double N_VDotProd_Serial(N_Vector x, N_Vector y)
{
  long int i, N;
  double sum, *xd, *yd;

  sum = ZERO;
  xd = yd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);

  for (i = 0; i < N; i++)
    sum += xd[i]*yd[i];
  
  return(sum);
}

double N_VMaxNorm_Serial(N_Vector x)
{
  long int i, N;
  double max, *xd;

  max = ZERO;
  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++) {
    if (fabs(xd[i]) > max) max = fabs(xd[i]);
  }

  return(max);
}

double N_VWrmsNorm_Serial(N_Vector x, N_Vector w)
{
  long int i, N;
  double sum, prodi, *xd, *wd;

  sum = ZERO;
  xd = wd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  wd = NV_DATA_S(w);

  for (i = 0; i < N; i++) {
    prodi = xd[i]*wd[i];
    sum += SQR(prodi);
  }

  return(RSqrt(sum/N));
}

double N_VWrmsNormMask_Serial(N_Vector x, N_Vector w, N_Vector id)
{
  long int i, N;
  double sum, prodi, *xd, *wd, *idd;

  sum = ZERO;
  xd = wd = idd = NULL;

  N  = NV_LENGTH_S(x);
  xd  = NV_DATA_S(x);
  wd  = NV_DATA_S(w);
  idd = NV_DATA_S(id);

  for (i = 0; i < N; i++) {
    if (idd[i] > ZERO) {
      prodi = xd[i]*wd[i];
      sum += SQR(prodi);
    }
  }

  return(RSqrt(sum / N));
}

double N_VMin_Serial(N_Vector x)
{
  long int i, N;
  double min, *xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  min = xd[0];

  for (i = 1; i < N; i++) {
    if (xd[i] < min) min = xd[i];
  }

  return(min);
}

void N_VCompare_Serial(double c, N_Vector x, N_Vector z)
{
  long int i, N;
  double *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) {
    zd[i] = (fabs(xd[i]) >= c) ? ONE : ZERO;
  }

  return;
}

int N_VInvTest_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  double *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++) {
    if (xd[i] == ZERO) return(FALSE);
    zd[i] = ONE/xd[i];
  }

  return(TRUE);
}

/*
 * -----------------------------------------------------------------
 * private functions
 * -----------------------------------------------------------------
 */

static void VCopy_Serial(N_Vector x, N_Vector z)
{
  long int i, N;
  double *xd, *zd;

  xd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]; 

  return;
}

static void VSum_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int i, N;
  double *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (i = 0; i < N; i++)
    zd[i] = xd[i]+yd[i];

  return;
}

static void VDiff_Serial(N_Vector x, N_Vector y, N_Vector z)
{
  long int N;
  double *xd, *yd, *zd;

  xd = yd = zd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  yd = NV_DATA_S(y);
  zd = NV_DATA_S(z);

  for (size_t i = 0; i < N; i++)
    zd[i] = xd[i]-yd[i];

  return;
}

static void VNeg_Serial(N_Vector x, N_Vector z) {
  long int N;
  double *xd, *zd;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);
  zd = NV_DATA_S(z);
  
  for (size_t i = 0; i < N; i++)
    zd[i] = -xd[i];

  return;
}



static void VScaleDiff_Serial(double c, N_Vector x, N_Vector y, N_Vector z)
{
    const long int N  = NV_LENGTH_S(x);
    const double *xd = NV_DATA_S(x);
    const double *yd = NV_DATA_S(y);
    double *zd = NV_DATA_S(z);

  for (size_t i = 0; i < N; i++)
    zd[i] = c*(xd[i]-yd[i]);

  return;
}

static void VLin1_Serial(double a, N_Vector x, N_Vector y, N_Vector z)
{
    const long int N  = NV_LENGTH_S(x);
    const double *xd = NV_DATA_S(x);
    const double *yd = NV_DATA_S(y);
    double *zd = NV_DATA_S(z);
    
    for (size_t i = 0; i < N; i++)
        zd[i] = (a*xd[i])+yd[i];
    
    return;
}

static void VLin2_Serial(double a, N_Vector x, N_Vector y, N_Vector z)
{
    const long int N  = NV_LENGTH_S(x);
    const double *xd = NV_DATA_S(x);
    const double *yd = NV_DATA_S(y);
    double * const zd = NV_DATA_S(z);
    
    for (size_t i = 0; i < N; i++)
        zd[i] = (a*xd[i])-yd[i];
    
    return;
}

static inline void Vaxpy_Serial(double a, N_Vector x, N_Vector y)
{
    const long int N  = NV_LENGTH_S(x);
    const double * const xd = NV_DATA_S(x);
    double * const yd = NV_DATA_S(y);
    
    for (size_t i = 0; i < N; i++)
        yd[i] += a*xd[i];
    
    return;
}

static void VScaleBy_Serial(double a, N_Vector x)
{
  long int i, N;
  double *xd;

  xd = NULL;

  N  = NV_LENGTH_S(x);
  xd = NV_DATA_S(x);

  for (i = 0; i < N; i++)
    xd[i] *= a;

  return;
}

/*
 *  TAMode.h
 *  TAMode
 *
 *  Created by Aaron Meyer on 2/15/15.
 *  Copyright (c) 2015 Aaron Meyer. All rights reserved.
 *
 */

#ifndef TAMode_XX
#define TAMode_XX

#include "sundials_nvector.h"

#define autocrineT 100000
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

struct TAMout {
    double *total;
    double *pY;
    double *surf;
};


int calcProfileMatlab(double *, double *, double *, int, double, int);
void *initState( N_Vector, struct rates *);

int calcSingleTAMdose (TAMout *, struct rates *, double, double *, unsigned int);

/* The classes below are exported */
extern "C" int calcProfileMat (double *, double *, int, double *, double *, double *, double);
extern "C" double calcLew (double *params);

#endif

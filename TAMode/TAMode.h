/*
 *  TAMode.h
 *  TAMode
 *
 *  Created by Aaron Meyer on 2/15/15.
 *  Copyright (c) 2015 Aaron Meyer. All rights reserved.
 *
 */

#ifndef TAMode_
#define TAMode_

/* The classes below are exported */
#pragma GCC visibility push(default)
double pyEntry(double *pIn);
int calcProfileMatlab(double *, double *, double *, int, double, int);
void pyEntryVec(double *pIn, double *pOut, int n);
double pyEntryNew(double *);

#pragma GCC visibility pop
#endif

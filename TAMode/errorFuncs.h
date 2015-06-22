//
//  errorFuncs.h
//  TAMode
//
//  Created by Aaron Meyer on 6/19/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#ifndef __TAMode__errorFuncs__
#define __TAMode__errorFuncs__

#include <stdio.h>

double errorFuncFix (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN);
double errorFuncOpt (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN, double *xx);

#endif /* defined(__TAMode__errorFuncs__) */

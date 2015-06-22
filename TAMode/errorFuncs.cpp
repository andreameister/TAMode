//
//  errorFuncs.cpp
//  TAMode
//
//  Created by Aaron Meyer on 6/19/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#include "errorFuncs.h"
#include <cmath>
#include "cobyla.h"
#include <string>
#include <exception>

using namespace std;

struct inData {
    size_t N;
    const double *fitt;      // Calculated model values.
    const double *pYmeas;    // pY measurement.
    const double *errorMeas; // Error for pY measurement.
};

static double errorFunc (double fitt, double pYmeas, double errorMeas) {
    return pow((((double) fitt) - pYmeas) / errorMeas, 2) / 2;
}

static double errorOpt(unsigned, const double *x, double *, void *data) {
    struct inData *dataS = (struct inData *) data;
    double xx = 0;
    
    for (int ii = 0; ii < dataS->N; ii++)
        xx += errorFunc(dataS->fitt[ii] * x[0], dataS->pYmeas[ii], dataS->errorMeas[ii]);
    
    return xx;
}

static double initialCondition (struct inData *dataS) {
    double meas = 0;
    double fit = 0;
    
    for (int ii = 0; ii < dataS->N; ii++) {
        meas += dataS->fitt[ii];
        fit += dataS->pYmeas[ii];
    }
    
    return fit / meas;
}

double errorFuncOpt (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN, double *xx) {
    struct inData dataS;
    dataS.fitt = fitt;
    dataS.pYmeas = pYmeas;
    dataS.errorMeas = errorMeas;
    dataS.N = inN;
    
    double ff = 0;
    xx[0] = initialCondition(&dataS);
    const double lower = xx[0]/2;
    const double upper = xx[0]*2;
    double dx = 3*xx[0]/8;
    double del = 1E-8;
    
    nlopt_stopping stop;
    stop.n = 0;
    stop.minf_max = 0.0;
    stop.ftol_rel = 1E-8;
    stop.ftol_abs = 0;
    stop.xtol_rel = 1E-8;
    stop.xtol_abs = &del;
    stop.nevals = 0;
    stop.maxeval = 1E9;
    stop.force_stop = 0;
    
    int flag = cobyla_minimize(1, errorOpt, &dataS, 0, NULL, 0, NULL, &lower, &upper, xx, &ff, &stop, &dx);
    
    if (flag < 0) throw runtime_error(string("Error during error optimization step."));
    
    return ff;
}

double errorFuncFix (const double *fitt, const double *pYmeas, const double *errorMeas, size_t inN) {
    double xx = 0;
    
    for (int ii = 0; ii < inN; ii++)
        xx += errorFunc(fitt[ii], pYmeas[ii], errorMeas[ii]);
    
    return xx;
}




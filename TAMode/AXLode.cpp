//
//  AXLode.cpp
//  TAMode
//
//  Created by Aaron Meyer on 6/19/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#include "AXLode.h"
#include "reactCode.h"
#include "sundials_nvector.h"
#include "CVodeHelpers.h"
#include "TAMode.h"
#include <string>
#include "errorFuncs.h"
#include <exception>

using namespace std;


static const double times[2] = {60, 240}; ///< Times of kinetic measurements.
static const double Gass[6] = {64, 16, 4, 1, 0.25, 0}; ///< Kinetic Gas6 doses.
static const double kTPS[5] = {0, 0.5, 1, 5, 10};


static const double pYk[5] = {4.1, 3.6, 5.3, 11.3, 11.6};
static const double pYkErr[5] = {1.0, 1.3, 1.1, 0.7, 0.4};


// Wrapping is outermost cell line, then Gas, then time
static const double pY[6][2] = { ///< pY measurements on short time scales.
    {10.75952427, 8.305264139},
    {7.390399159, 7.056438019},
    {7.144036441, 7.680851079},
    {4.570826833, 8.184089069},
    {6.107714557, 7.204021903},
    {7.535575387, 7.535575387}};

static const double pYerror[6][2] = { ///< Error for short time scale pY measurements.
    {1.74431775,  2.100723242},
    {1.267611,    1.260108508},
    {0.898008437, 1.680415875},
    {1.521844479, 0.871927763},
    {0.932623012, 0.563182873},
    {0.812417951, 0.812417951}};

static const double tot[6][2] = {
    {3443.11, 3219.69},
    {3143.41, 3353.82},
    {3018.88, 3611.82},
    {2608.88, 3448.21},
    {2690.24, 3168.14},
    {2672.00, 2672.00}}; // A549

static const double totError[6][2] = {
    {174.38, 132.10},
    {189.03, 129.93},
    {245.75, 225.42},
    {154.89, 203.72},
    {128.72, 187.34},
    {82.62, 82.62}}; // A549

static const double surf[6][2] = {
    {0.206, 0.239},
    {0.274, 0.316},
    {0.281, 0.251},
    {0.220, 0.302},
    {0.256, 0.281},
    {0.257, 0.337}}; // A549

static const double surfError[6][2] = {
    {0.043, 0.015},
    {0.047, 0.037},
    {0.032, 0.030},
    {0.025, 0.036},
    {0.044, 0.035},
    {0.030, 0.023}}; // A549


static void calcKinetic (double *outData, double *totData, double *surfData, double *earlyPY, struct rates *params) {
    N_Vector init_state = N_VNew_Serial(Nspecies);
    double t;
    
    void *cvode_mem = initState(init_state, params);
    // Initialize state based on autocrine ligand
    
    if (cvode_mem == NULL) {
        N_VDestroy_Serial(init_state);
        throw runtime_error(string("Error during solver threads."));
        return;
    }
    
    //
    // This part will calculate the dose response
    //
    struct rates paramTwo = *params;
    N_Vector state = N_VClone(init_state);
    // Initialize state based on autocrine ligand
    
    for (size_t stimuli = 0; stimuli < 6; stimuli++) {
        for (int xx = 0; xx < Nspecies; xx++) Ith(state,xx) = Ith(init_state,xx);
        
        paramTwo.gasCur = paramTwo.autocrine + Gass[stimuli];
        CVodeSetUserData(cvode_mem, &paramTwo);
        
        t = 0;
        
        solverReset(cvode_mem, state);
        
        /* In loop, call CVode, print results, and test for error.
         Break out of loop when NOUT preset output times have been reached.  */
        
        for (unsigned int ii = 0; ii < NELEMS(times); ii++) {
            int flag = CVode(cvode_mem, times[ii], state, &t, CV_NORMAL);
            
            if (flag < 0) {
                N_VDestroy_Serial(state);
                CVodeFree(&cvode_mem);
                N_VDestroy_Serial(init_state);
                throw runtime_error(string("Error during solver threads."));
                return;
            }
            
            outData[stimuli*NELEMS(times) + ii] = pYCalc(state,params,AXL);
            totData[stimuli*NELEMS(times) + ii] = totCalc(state,params);
            surfData[stimuli*NELEMS(times) + ii] = surfCalc(state, AXL);
        }
    }
    
    //
    // This part calculates the kinetic response
    //
    earlyPY[0] = pYCalc(init_state,params,AXL);
    // Initialize state based on autocrine ligand
    
    /* We've got the initial state, so now run through the kinetic data */
    for (int xx = 0; xx < Nspecies; xx++) Ith(state,xx) = Ith(init_state,xx);
    
    paramTwo.gasCur = paramTwo.autocrine + 1.25;
    
    CVodeSetUserData(cvode_mem, &paramTwo);
    
    t = 0;
    
    solverReset(cvode_mem, state);
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    
    for (unsigned int ii = 1; ii < NELEMS(kTPS); ii++) {
        int flag = CVode(cvode_mem, kTPS[ii], state, &t, CV_NORMAL);
        
        if (flag < 0) {
            N_VDestroy_Serial(state);
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(init_state);
            throw runtime_error(string("Error during solver threads."));
            return;
        }
        
        earlyPY[ii] = pYCalc(state,params,AXL);
    }
    
    N_VDestroy_Serial(state);
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(init_state);
}



double calcErrorA549 (struct rates inP) {
    double outData[NELEMS(Gass)*NELEMS(times)];
    double totData[NELEMS(Gass)*NELEMS(times)];
    double surfData[NELEMS(Gass)*NELEMS(times)];
    double earlyPY[NELEMS(kTPS)];
    
    double fitParam;
    
    double error = 0;
    
    try {
        calcKinetic(outData, totData, surfData, earlyPY, &inP);
        
        error += errorFuncOpt (outData, pY[0], pYerror[0], NELEMS(outData), &fitParam);
        error += errorFuncFix (totData, tot[0], totError[0], NELEMS(totData));
        error += errorFuncOpt (surfData, surf[0], surfError[0], NELEMS(surfData), &fitParam);
        error += errorFuncOpt (earlyPY, pYk, pYkErr, NELEMS(earlyPY), &fitParam);
    } catch (runtime_error &e) {
        errorLogger(&e);
        error = 1E8;
    }
    
    return error;
}

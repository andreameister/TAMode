#include <iostream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include "TAMode.h"
#include "sundials_nvector.h"
#include "cvode.h"
#include <vector>
#include <string>
#include <exception>
#include <cmath>
#include "cvode_impl.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "CVodeHelpers.h"
#include "reactCode.h"

using namespace std;


// Calculate the initial state by waiting a long time with autocrine Gas
void *initState( N_Vector init, struct rates *params) {
    double t;
    
    for (int ii = 0; ii < Nspecies ; ii++) Ith(init,ii) = 0;
    
    params->gasCur = params->autocrine;
    
    void *cvode_mem = solver_setup (init, params, AXL_react);
    if (cvode_mem == NULL) return NULL;
    
    int flag = CVode(cvode_mem, autocrineT, init, &t, CV_NORMAL);
    if (flag < 0) {
        CVodeFree(&cvode_mem);
        return NULL;
    }
    
    /* Free integrator memory */
    return cvode_mem;
}

/// Calculate phosphorylation at time points measured
static void calcProfileSet (TAMout *outData, double *tps, struct rates *params, int nTps, double GasStim) {
    N_Vector state = N_VNew_Serial(Nspecies);
    struct rates paramTwo = *params;
    
    double t; ///< Time position of the solver.
    
    void *cvode_mem = NULL;
    int flag;
    
    // Initialize state based on autocrine ligand
    try {
        cvode_mem = initState(state, &paramTwo);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        throw;
    }
    
    
    /* We've got the initial state, so now run through the kinetic data */
    paramTwo.gasCur = paramTwo.gasCur + GasStim;
    CVodeSetUserData(cvode_mem, &paramTwo);
    t = 0;
    
    try {
        solverReset(cvode_mem, state);
    } catch (exception &e) {
        N_VDestroy_Serial(state);
        CVodeFree(&cvode_mem);
        throw;
    }
    
    /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
    size_t ii = 0;
    
    if (tps[0] == 0) {
        outData->total[0] = receptorTotCalc(state, params, AXL);
        outData->total[1] = receptorTotCalc(state, params, Mer);
        outData->total[2] = receptorTotCalc(state, params, Tyro);
        
        outData->pY[0] = pYCalc(state, params, AXL);
        outData->pY[1] = pYCalc(state, params, Mer);
        outData->pY[2] = pYCalc(state, params, Tyro);
        
        outData->surf[0] = surfCalc(state, AXL);
        outData->surf[1] = surfCalc(state, Mer);
        outData->surf[2] = surfCalc(state, Tyro);
        
        ii = 1;
    }
    
    for (; ii < (size_t) abs(nTps); ii++) {
        flag = CVode(cvode_mem, tps[ii], state, &t, CV_NORMAL);
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(state);
            throw runtime_error(string("Error at CVode Time Course."));
        }
        
        outData->total[3*ii] = receptorTotCalc(state, params, AXL);
        outData->total[3*ii + 1] = receptorTotCalc(state, params, Mer);
        outData->total[3*ii + 2] = receptorTotCalc(state, params, Tyro);
        
        outData->pY[3*ii] = pYCalc(state, params, AXL);
        outData->pY[3*ii + 1] = pYCalc(state, params, Mer);
        outData->pY[3*ii + 2] = pYCalc(state, params, Tyro);
        
        outData->surf[3*ii] = surfCalc(state, AXL);
        outData->surf[3*ii + 1] = surfCalc(state, Mer);
        outData->surf[3*ii + 2] = surfCalc(state, Tyro);
    }
    
    /* Free integrator memory */
    CVodeFree(&cvode_mem);
    N_VDestroy_Serial(state);
}




/// Single receptor dose response
int calcSingleTAMdose (TAMout *outData, struct rates *params, double tp, double *doses, unsigned int nDoses) {
    struct rates inP;
    N_Vector init = N_VNew_Serial(Nspecies);
    N_Vector state = N_VNew_Serial(Nspecies);
    
    void *cvode_mem = NULL;
    int flag;
    double t; ///< Time position of the solver.
    
    for (size_t recp = 0; recp < 3; recp++) {
        inP = *params;
        
        for (size_t jj = 0; jj < 3; jj++) {
            if (jj != recp) inP.TAMs[jj].expression = 0;
        }
        
        try {
            cvode_mem = initState(init, &inP);
        } catch (exception &e) {
            CVodeFree(&cvode_mem);
            N_VDestroy_Serial(init);
            N_VDestroy_Serial(state);
            return -1;
        }
        
    
        for (size_t ii = 0; ii < nDoses; ii++) {
            inP.gasCur = inP.autocrine + doses[ii];
            CVodeSetUserData(cvode_mem, &inP);
            t = 0;
            
            for (int xx = 0; xx < Nspecies ; xx++) Ith(state,xx) = Ith(init,xx);
            
            
            try {
                solverReset(cvode_mem, state);
            } catch (exception &e) {
                N_VDestroy_Serial(state);
                N_VDestroy_Serial(init);
                CVodeFree(&cvode_mem);
                return -1;
            }
            
            if (tp > 0) {
                flag = CVode(cvode_mem, tp, state, &t, CV_NORMAL);
                if (flag < 0) {
                    CVodeFree(&cvode_mem);
                    N_VDestroy_Serial(state);
                    N_VDestroy_Serial(init);
                    return -1;
                }
            }
            
            if (recp == 0) {
                outData->pY[recp*nDoses + ii] = pYCalc(state, params, AXL);
                outData->surf[recp*nDoses + ii] = surfCalc(state, AXL);
                outData->total[recp*nDoses + ii] = receptorTotCalc(state, params, AXL);
            } else if (recp == 1) {
                outData->pY[recp*nDoses + ii] = pYCalc(state, params, Mer);
                outData->surf[recp*nDoses + ii] = surfCalc(state, Mer);
                outData->total[recp*nDoses + ii] = receptorTotCalc(state, params, Mer);
            } else {
                outData->pY[recp*nDoses + ii] = pYCalc(state, params, Tyro);
                outData->surf[recp*nDoses + ii] = surfCalc(state, Tyro);
                outData->total[recp*nDoses + ii] = receptorTotCalc(state, params, Tyro);
            }
        }
        
        CVodeFree(&cvode_mem);
    }
    
    N_VDestroy_Serial(state);
    N_VDestroy_Serial(init);
    return 0;
}

static double LewError (double *in1, double *in2, unsigned int len) {
    vector<double> temp1(in1, in1+len);
    vector<double> temp2(in2, in2+len);
    
    const double mag1 = accumulate(temp1.begin(), temp1.end(), (double) 0.0);
    const double mag2 = accumulate(temp2.begin(), temp2.end(), (double) 0.0);
    
    for (unsigned int ii = 0; ii < len; ii++) {
        temp1[ii] = fabs((temp1[ii] / mag1) - (temp2[ii] / mag2));
    }
    
    return accumulate(temp1.begin(), temp1.end(), (double) 0.0);
}

extern "C" double calcLew (double *inP) {
    struct rates params = Param(inP);
    
    TAMout outData;
    
    outData.pY = (double *) malloc(sizeof(double)*3*6);
    outData.total = (double *) malloc(sizeof(double)*3*6);
    outData.surf = (double *) malloc(sizeof(double)*3*6);
    
    double doses[] = {0, 25, 50, 100, 200, 300};
    
    double TyroPY[] = {1, 1.5, 2.0, 5, 10, 10};
    double AxlPY[] = {1, 1, 1, 1, 1, 1};
    double MerPY[] = {0, 1, 3, 5, 10, 10};
    
    int flag = calcSingleTAMdose (&outData, &params, 10, doses, NELEMS(doses));
    
    if (flag < 0) {
        free(outData.pY);
        free(outData.total);
        free(outData.surf);
        return -1;
    }
    
    double error = 0;
    
    error += LewError(AxlPY, &outData.pY[0], NELEMS(doses));
    error += LewError(MerPY, &outData.pY[NELEMS(doses)], NELEMS(doses));
    error += LewError(TyroPY, &outData.pY[2*NELEMS(doses)], NELEMS(doses));
    
    free(outData.pY);
    free(outData.total);
    free(outData.surf);
    
    return error;
}

extern "C" int calcProfileMat (double *paramIn, double *tps, int nTps, double *total, double *pY, double *surf, double GasStim) {
    TAMout data;
    data.total = total;
    data.pY = pY;
    data.surf = surf;
    
    rates r = Param(paramIn);
    
    try {
        calcProfileSet (&data, tps, &r, nTps, GasStim);
        return 0;
    } catch (exception &e) {
        return -1;
    }
}

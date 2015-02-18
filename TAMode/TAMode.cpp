#include <iostream>
#include <sstream>
#include "TAMode.h"
#include "CVodeHelpers.h"


using namespace std;

//double pyEntry(double *pIn) {
//    return calcError(Param(pIn));
//}
//
//void pyEntryVec(double *pIn, double *pOut, int n) {
//    for (int ii = 0; ii < n; ii++) {
//        if (pIn[ii*11 + 0] < pIn[ii*11 + 2]) {
//            pOut[ii] = 1E6;
//            continue;
//        }
//        
//        pOut[ii] = calcError(Param(&pIn[ii*11]));
//    }
//}
//
//int calcProfileMatlab(double *dataPtr, double *params, double *tps, int nTps, double GasStim, int frac) {
//    struct rates pInS = Param(params);
//    
//    try {
//        calcProfileSet (dataPtr, tps, &pInS, nTps, GasStim, frac);
//    } catch (std::exception &e) {
//        errorLogger(&e);
//        return 1;
//    }
//    
//    return 0;
//}



//
//  ModelRunning.cpp
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#include "sundials_nvector.h"
#include "cvode.h"
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <cmath>
#include "cvode_impl.h"
#include "CVodeHelpers.h"

#define NVp N_VGetArrayPointer_Serial



using namespace std;

//const double fgMgConv = 135.2;

void trafFunc(double *dextR, double *dintR, double intRate, struct rates *tr, double extR, double intR) {
    *dextR += -extR*intRate + tr->kRec*(1-tr->fElse)*intR*tr->internalFrac; // Endocytosis, recycling
    *dintR += extR*intRate/tr->internalFrac - tr->kRec*(1-tr->fElse)*intR - tr->kDeg*tr->fElse*intR; // Endocytosis, recycling, degradation
}

void heteroTAM (double *Rone, double *Rtwo, double *dRone, double *dRtwo, hetRates *hetR, double *hetDim, double *dhetDim, struct rates *tr) {
    const double dR7 = hetR->xFwd7*Rone[2]*Rtwo[2] - hetR->xRev7*hetDim[2];
    const double dR8 = hetR->xFwd8*Rone[0]*Rtwo[3] - hetR->xRev8*hetDim[2];
    const double dR9 = hetR->xFwd9*Rone[1]*Rtwo[1] - hetR->xRev9*hetDim[2];
    const double dR10 = hetR->xFwd10*Rone[3]*Rtwo[0] - hetR->xRev10*hetDim[2];
    const double dR11 = hetR->xFwd11*Rone[2]*Rtwo[0] - hetR->xRev11*hetDim[1];
    const double dR12 = hetR->xFwd12*Rone[0]*Rtwo[1] - hetR->xRev12*hetDim[1];
    const double dR13 = hetR->xFwd13*Rone[1]*Rtwo[0] - hetR->xRev13*hetDim[0];
    const double dR14 = hetR->xFwd14*Rone[0]*Rtwo[2] - hetR->xRev14*hetDim[0];
    
    const double dR7i = hetR->xFwd7*Rone[8]*Rtwo[8] - hetR->xRev7*hetDim[5];
    const double dR8i = hetR->xFwd8*Rone[6]*Rtwo[9] - hetR->xRev8*hetDim[5];
    const double dR9i = hetR->xFwd9*Rone[7]*Rtwo[7] - hetR->xRev9*hetDim[5];
    const double dR10i = hetR->xFwd10*Rone[9]*Rtwo[6] - hetR->xRev10*hetDim[5];
    const double dR11i = hetR->xFwd11*Rone[8]*Rtwo[6] - hetR->xRev11*hetDim[4];
    const double dR12i = hetR->xFwd12*Rone[8]*Rtwo[7] - hetR->xRev12*hetDim[4];
    const double dR13i = hetR->xFwd13*Rone[7]*Rtwo[6] - hetR->xRev13*hetDim[3];
    const double dR14i = hetR->xFwd14*Rone[6]*Rtwo[8] - hetR->xRev14*hetDim[3];
    
    dRone[0] += -dR8 - dR12 - dR14;
    dRone[1] += -dR9 - dR13;
    dRone[2] += -dR7 - dR11;
    dRone[3] += -dR10;
    dRone[6] += -dR8i - dR12i - dR14i;
    dRone[7] += -dR9i - dR13i;
    dRone[8] += -dR7i - dR11i;
    dRone[9] += -dR10i;
    
    
    dRtwo[0] += -dR10 - dR11 - dR13;
    dRtwo[1] += -dR9 - dR12;
    dRtwo[2] += -dR7 - dR14;
    dRtwo[3] += -dR8;
    dRtwo[6] += -dR10i - dR11i - dR13i;
    dRtwo[7] += -dR9i - dR12i;
    dRtwo[8] += -dR7i - dR14i;
    dRtwo[9] += -dR8i;
    
    dhetDim[0] += dR13 + dR14; // AMD1
    dhetDim[1] += dR11 + dR12; // MAD1
    dhetDim[2] += dR7 + dR8 + dR9 + dR10; // AMD2
    
    dhetDim[3] += dR13i + dR14i;
    dhetDim[4] += dR11i + dR12i;
    dhetDim[5] += dR7i + dR8i + dR9i + dR10i;
    
    trafFunc(&dhetDim[0], &dhetDim[3], tr->internalize + (tr->pYinternalize*((double) tr->pD1 == 1)), tr, hetDim[0], hetDim[3]);
    trafFunc(&dhetDim[1], &dhetDim[4], tr->internalize + (tr->pYinternalize*((double) tr->pD1 == 1)), tr, hetDim[4], hetDim[10]);
    trafFunc(&dhetDim[2], &dhetDim[5], tr->internalize + tr->pYinternalize, tr, hetDim[4], hetDim[10]);
}



void TAM_react(double *R, double *Li, double *dR, double *dLi, struct TAMrates *r, struct rates *tr) {
    const double dR1 = r->Binding1 * R[0] * tr->gasCur - r->Unbinding1 * R[1];
    const double dR2 = r->Binding2 * R[0] * tr->gasCur - r->Unbinding2 * R[2];
    const double dR3 = r->Binding2 * R[1] * tr->gasCur - r->Unbinding2 * R[3];
    const double dR4 = r->Binding1 * R[2] * tr->gasCur - r->Unbinding1 * R[3];
    const double dR5 = r->xFwd1 * R[0] * R[1] - r->xRev1 * R[4];
    const double dR6 = r->xFwd2 * R[0] * R[2] - r->xRev2 * R[4];
    const double dR7 = r->xFwd3 * R[0] * R[3] - r->xRev3 * R[5];
    const double dR8 = r->xFwd4 * R[1] * R[1] - r->xRev4 * R[5];
    const double dR9 = r->xFwd5 * R[2] * R[2] - r->xRev5 * R[5];
    const double dR11 = r->xFwd6 * tr->gasCur * R[4] - r->xRev6 * R[5];
    
    const double dR32 = r->Binding1 * R[6] * R[12] / tr->internalV - r->Unbinding1 * R[7];
    const double dR33 = r->Binding2 * R[6] * R[12] / tr->internalV - r->Unbinding2 * R[8];
    const double dR34 = r->Binding2 * R[7] * R[12] / tr->internalV - r->Unbinding2 * R[9];
    const double dR35 = r->Binding1 * R[8] * R[12] / tr->internalV - r->Unbinding1 * R[9];
    const double dR36 = r->xFwd1 * R[6] * R[7] - r->xRev1 * R[10];
    const double dR37 = r->xFwd2 * R[6] * R[8] - r->xRev2 * R[10];
    const double dR38 = r->xFwd3 * R[6] * R[9] - r->xRev3 * R[11];
    const double dR39 = r->xFwd4 * R[7] * R[7] - r->xRev4 * R[11]; // Checked
    const double dR40 = r->xFwd5 * R[8] * R[8] - r->xRev5 * R[11]; // Checked
    const double dR41 = r->xFwd6 * (*Li) * R[10] / tr->internalV - r->xRev6 * R[11]; // Checked
    
    dR[0] += - dR7 - dR6 - dR5 - dR1 - dR2 + r->expression; // AXL
    dR[1] += -2*(dR8) - dR5 + dR1 - dR3; // AXLgas1
    dR[2] += -2*(dR9) - dR6 + dR2 - dR4; // AXLgas2
    dR[3] += -dR7 + dR3 + dR4; // AXLgas12
    dR[4] += -dR11 + dR6 + dR5; // AXLdimer1
    dR[5] += dR11 + dR9 + dR8 + dR7; // AXLdimer2
    
    dR[6] += - dR38 - dR37 - dR36 - dR32 - dR33; // AXLi
    dR[7] += -2*(dR39) - dR36 + dR32 - dR34; // AXLgas1i
    dR[8] += -2*(dR40) - dR37 + dR33 - dR35; // AXLgas2i
    dR[9] += -dR38 + dR34 + dR35; // AXLgas12i
    dR[10] += -dR41 + dR37 + dR36; // AXLdimer1i
    dR[11] += dR41 + dR40 + dR39 + dR38; // AXLdimer2i
    
    *dLi += -dR41 - dR32 - dR33 - dR34 - dR35;
    
    for (int ii = 0; ii < 4; ii++) {
        trafFunc(&dR[ii], &dR[ii+6], tr->internalize, tr, R[ii], R[ii+6]);
    }
    
    trafFunc(&dR[4], &dR[10], tr->internalize + (tr->pYinternalize*((double) tr->pD1 == 1)), tr, R[4], R[10]);
    trafFunc(&dR[5], &dR[11], tr->internalize + tr->pYinternalize, tr, R[5], R[11]);
}

int AXL_react(double t, N_Vector xIn, N_Vector dxdtIn, void *user_data) {
    double* x_d = NV_DATA_S(xIn);
    double* dxdt_d = NV_DATA_S(dxdtIn);
    struct rates *r = (struct rates *) user_data;
    
    // 0 AXL   // 1 A1    // 2 A2
    // 3 A12    // 4 D1    // 5 D2    // 6 AXLi
    // 7 A1i    // 8 A2i   // 9 A12i // 10 D1i    // 11 D2i   // 12 Gasi
    
    TAM_react(&x_d[0], &x_d[12], &dxdt_d[0], &dxdt_d[12], &r->TAMs[0], r);
    TAM_react(&x_d[13], &x_d[12], &dxdt_d[13], &dxdt_d[12], &r->TAMs[1], r);
    TAM_react(&x_d[25], &x_d[12], &dxdt_d[25], &dxdt_d[12], &r->TAMs[2], r);
    
    // AM, MT, AT
    heteroTAM (&x_d[0], &x_d[13], &dxdt_d[0], &dxdt_d[13], &r->hetR[0], &x_d[37], &dxdt_d[37], r);
    heteroTAM (&x_d[13], &x_d[25], &dxdt_d[13], &dxdt_d[25], &r->hetR[1], &x_d[40], &dxdt_d[40], r);
    heteroTAM (&x_d[0], &x_d[25], &dxdt_d[0], &dxdt_d[25], &r->hetR[2], &x_d[43], &dxdt_d[43], r);
    
    
    dxdt_d[12] = -r->kDeg*x_d[12];
    
    
    return 0;
}

struct rates Param(double *params) {
    struct rates out;
    
    if (min_element(params,params+9) < 0) {
        throw invalid_argument(string("Parameter outside the physical range."));
    }
    
    out.internalize = 0.03;
    out.pYinternalize = params[0];
    out.fElse = params[1];
    out.autocrine = params[2];
    out.pD1 = (int) params[3];
    double xFwd = params[4];
    
    out.hetR[0].xRev7 = params[5];
    out.hetR[1].xRev7 = params[6];
    out.hetR[2].xRev7 = params[7];
    
    out.fD2 = 1;
    out.kRec = 5.8E-2;
    out.kDeg = 2.2E-3;
    out.internalFrac = 0.5;
    out.internalV = 623;
    
    out.TAMs[0].Binding1 = 1.2;
    out.TAMs[0].Binding2 = 0.06;
    out.TAMs[0].Unbinding1 = 0.042;
    out.TAMs[0].Unbinding2 = params[8];
    out.TAMs[0].xRev4 = params[9];
    out.TAMs[0].xFwd6 = params[10];
    out.TAMs[0].expression = params[11];
    
    out.TAMs[1].Binding1 = 0.06;
    out.TAMs[1].Binding2 = 0.06;
    out.TAMs[1].Unbinding1 = params[12];
    out.TAMs[1].Unbinding2 = params[13];
    out.TAMs[1].xRev4 = params[14];
    out.TAMs[1].xFwd6 = params[15];
    out.TAMs[1].expression = params[16];
    
    out.TAMs[2].Binding1 = 0.06;
    out.TAMs[2].Binding2 = 0.06;
    out.TAMs[2].Unbinding1 = params[17];
    out.TAMs[2].Unbinding2 = params[18];
    out.TAMs[2].xRev4 = params[19];
    out.TAMs[2].xFwd6 = params[20];
    out.TAMs[2].expression = params[21];
    
    // Detailed balance
    for (size_t ii = 0; ii < 3; ii++) {
        const double KD1 = out.TAMs[ii].Unbinding1/out.TAMs[ii].Binding1;
        const double KD2 = out.TAMs[ii].Unbinding2/out.TAMs[ii].Binding2;
        
        out.TAMs[ii].xFwd1 = xFwd;
        out.TAMs[ii].xFwd2 = xFwd;
        out.TAMs[ii].xFwd3 = xFwd;
        out.TAMs[ii].xFwd4 = xFwd;
        out.TAMs[ii].xFwd5 = xFwd;
        
        out.TAMs[ii].xRev1 = out.TAMs[ii].Unbinding2; // Assuming off rate of one ligand dimer is represented by ligand off rate.
        out.TAMs[ii].xRev2 = KD1*out.TAMs[ii].xFwd2*out.TAMs[ii].xRev1/KD2/out.TAMs[ii].xFwd1;
        out.TAMs[ii].xRev5 = out.TAMs[ii].xFwd5*out.TAMs[ii].xRev4/out.TAMs[ii].xFwd4*KD1*KD1/KD2/KD2;
        out.TAMs[ii].xRev6 = out.TAMs[ii].xFwd1*out.TAMs[ii].xFwd6*out.TAMs[ii].xRev5*KD2*KD2/out.TAMs[ii].xRev1/KD1/out.TAMs[ii].xFwd5;
        out.TAMs[ii].xRev3 = KD1*out.TAMs[ii].xFwd3*out.TAMs[ii].xRev4/out.TAMs[ii].xFwd4/KD2;
    }
    
    for (size_t ii = 0; ii < 3; ii++) {
        const double KD11 = out.TAMs[hetRi[ii][0]].Unbinding1 / out.TAMs[hetRi[ii][0]].Binding1;
        const double KD12 = out.TAMs[hetRi[ii][1]].Unbinding1 / out.TAMs[hetRi[ii][1]].Binding1;
        const double KD21 = out.TAMs[hetRi[ii][0]].Unbinding2 / out.TAMs[hetRi[ii][0]].Binding2;
        const double KD22 = out.TAMs[hetRi[ii][1]].Unbinding2 / out.TAMs[hetRi[ii][1]].Binding2;
        
        
        out.hetR[ii].xFwd7 = xFwd;
        out.hetR[ii].xFwd8 = xFwd;
        out.hetR[ii].xFwd9 = xFwd;
        out.hetR[ii].xFwd10 = xFwd;
        out.hetR[ii].xFwd11 = xFwd;
        out.hetR[ii].xFwd12 = xFwd;
        out.hetR[ii].xFwd13 = xFwd;
        out.hetR[ii].xFwd14 = xFwd;
        
        out.hetR[ii].xRev8 = out.hetR[ii].xRev7*out.hetR[ii].xFwd8*KD21/KD12/out.hetR[ii].xFwd7;
        out.hetR[ii].xRev9 = KD22*KD21/KD11/KD12*out.hetR[ii].xFwd9*out.hetR[ii].xRev7/out.hetR[ii].xFwd7;
        out.hetR[ii].xRev10 = out.hetR[ii].xRev7*out.hetR[ii].xFwd10*KD22/KD11/out.hetR[ii].xFwd7;
        
        // TODO: Detailed balance for these. At the moment assuming unbinding of ligand is equivalent
        out.hetR[ii].xRev11 = out.TAMs[hetRi[ii][1]].Unbinding1;
        out.hetR[ii].xRev12 = out.TAMs[hetRi[ii][0]].Unbinding2;
        out.hetR[ii].xRev13 = out.TAMs[hetRi[ii][1]].Unbinding2;
        out.hetR[ii].xRev14 = out.TAMs[hetRi[ii][0]].Unbinding1;
        
        // still need ligand binding to the one ligand dimer
    }
    
    return out;
}




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





/// END REACTION CODE












//
//  ModelRunning.h
//  UniformOptimization
//
//  Created by Aaron Meyer on 3/13/14.
//  Copyright (c) 2014 Aaron Meyer. All rights reserved.
//

#ifndef __UniformOptimization__ModelRunning__
#define __UniformOptimization__ModelRunning__

#include "sundials_nvector.h"
#include "cvode_impl.h"
#include "cvode.h"
#include "sundials_dense.h"

#define autocrineT 10000
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */
#define NELEMS(x)  (sizeof(x) / sizeof(x[0]))

#define Nspecies 13
#define Ith(v,i)    NV_Ith_S(v,i)       /* Ith numbers components 1..NEQ */


struct TAMrates {
    double Binding1;   ///< Forward binding rate for Ig1
    double Binding2;   ///< Forward binding rate for Ig2
    double Unbinding1; ///< Reverse binding rate for Ig1
    double Unbinding2; ///< Reverse binding rate for Ig2
    double xFwd1;      ///< Reaction 1 forward rate.
    double xRev1;      ///< Reaction 1 reverse rate.
    double xFwd3;      ///< Reaction 3 forward rate.
    double xRev3;      ///< Reaction 3 reverse rate.
    double expression; ///< AXL expression rate.
    
    double xRev5;
    double xRev4;
    double xRev2;
    double xFwd2;
    double xFwd4;
    double xFwd5;
    double xFwd6;
    double xRev6;
};

struct hetRates {
    double xFwd7;
    double xRev7;
    double xFwd8;
    double xRev8;
    double xFwd9;
    double xRev9;
    double xFwd10;
    double xRev10;
    double xFwd11;
    double xRev11;
    double xFwd12;
    double xRev12;
    double xFwd13;
    double xRev13;
    double xFwd14;
    double xRev14;
};

struct rates {
    TAMrates TAMs[3]; // AXL, MerTK, Tyro3

    double kRec;       ///< Recycling rate.
    double kDeg;       ///< Degradation rate.
    double fElse;      ///< Recycling fraction for non-D2 species.
    double fD2;        ///< Recycling fraction for D2.
    double internalize;///< Non-pY species internalization rate.
    double pYinternalize;///< pY species internalization rate.
    double internalFrac;
    double internalV;
    double autocrine;
    double gasCur;
    int pD1;
    
    hetRates hetR[3]; // AM, MT, AT
    
    size_t hetRi[3][2] = {{1, 2}, {2, 3}, {1, 3}};
};



int AXL_react(double, N_Vector, N_Vector, void *);
struct rates Param(double*);

#endif /* defined(__UniformOptimization__ModelRunning__) */

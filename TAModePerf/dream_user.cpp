//
//  main.cpp
//  TAModePerf
//
//  Created by Aaron Meyer on 3/11/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "cvode.h"
#include <time.h>
#include "TAMode.h"
#include "dream_user.h"
#include <algorithm>
#include "pdflib.hpp"
#include <atomic>

using namespace std;

atomic<double> maxVal;


const double lower[7] = {-10, -5, 0, -5, -5, -5, -5};
const double upper[7] = {5,    5, 5,  5,  5,  5,  5};

double sample_likelihood ( int, double zp[] ) {
    double zpL[NELEMS(lower)];
    
    for (size_t ii = 0; ii < NELEMS(lower); ii++) {
        zpL[ii] = pow(10,zp[ii]);
    }
    
    const double outter = -calcLew (zpL);
    
    if (outter > (maxVal + 0.1)) {
        cout << outter << endl;
        maxVal = outter;
    }
    
    return outter;
}

//  Parameters:
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, double PRIOR_SAMPLE[PAR_NUM], the sample from the distribution.
double *prior_sample ( int par_num ) {
    double *zp;
    
    zp = (double *) malloc((unsigned long) par_num * sizeof(double));
    
    for (size_t i = 0; i < par_num; i++ ) {
        zp[i] = r8_uniform_sample (lower[i], upper[i]);
    }
    
    return zp;
}

//  Parameters:
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Input, double ZP[PAR_NUM], the argument of the density
//    function.
//
//    Output, real PRIOR_DENSITY, the value of the prior density function.
double prior_density ( int par_num, double zp[] ) {
    double value;
    
    value = 1.0;
    
    for (size_t i = 0; i < par_num; i++ ) {
        value = value * r8_uniform_pdf ( lower[i], upper[i], zp[i] );
    }
    
    return value;
}

void problem_size ( int &chain_num, int &cr_num, int &gen_num, int &pair_num,
                   int &par_num ) {
//  Parameters:
//
//    Output, int &CHAIN_NUM, the total number of chains.
//    3 <= CHAIN_NUM.
//
//    Output, int &CR_NUM, the total number of CR values.
//    1 <= CR_NUM.
//
//    Output, int &GEN_NUM, the total number of generations.
//    2 <= GEN_NUM.
//
//    Output, int &PAIR_NUM, the number of pairs of
//    crossover chains.
//    0 <= PAIR_NUM.
//
//    Output, int &PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
    chain_num = 8;
    cr_num = 5;
    gen_num = 1E5;
    pair_num = 2;
    par_num = NELEMS(lower);
    
    return;
}



void problem_value ( string *chain_filename, string *gr_filename,
                    double &gr_threshold, int &jumpstep, double limits[], int par_num,
                    int &printstep, string *restart_read_filename,
                    string *restart_write_filename )
{
maxVal = -10000;
//  Parameters:
//
//    Output, string CHAIN_FILENAME, the "base" filename
//    to be used for the chain files.  If this is ""
//    then the chain files will not be written.  This name should
//    include a string of 0's which will be replaced by the chain
//    indices.  For example, "chain000.txt" would work as long as the
//    number of chains was 1000 or less.
//
//    Output, string *GR_FILENAME, the name of the file
//    in which values of the Gelman-Rubin statistic will be recorded,
//    or "" if this file is not to be written.
//
//    Output, double &GR_THRESHOLD, the convergence tolerance for
//    the Gelman-Rubin statistic.
//
//    Output, int &JUMPSTEP, forces a "long jump" every
//    JUMPSTEP generations.
//
//    Output, double LIMITS[2*PAR_NUM], lower and upper bounds
//    for each parameter.
//
//    Input, int PAR_NUM, the total number of parameters.
//    1 <= PAR_NUM.
//
//    Output, int &PRINTSTEP, the interval between generations on
//    which the Gelman-Rubin statistic will be computed and written to a file.
//
//    Output, string *RESTART_READ_FILENAME, the name of the file
//    containing restart information.  If this calculation is not a restart,
//    then this should be "".
//
//    Output, string *RESTART_WRITE_FILENAME, the name of the file
//    to be written, containing restart information.  If a restart file is not
//    to be written, this should be "".
//

    *chain_filename = "problem0_chain00.txt";
    *gr_filename = "problem0_gr.txt";
    gr_threshold = 1.1;
    jumpstep = 5;
    for (size_t j = 0; j < par_num; j++ ) {
        limits[0+j*2] = lower[j];
        limits[1+j*2] = upper[j];
    }
    printstep = 30;
    *restart_read_filename = "";
    *restart_write_filename = "problem0_restart.txt";
    
    return;
}

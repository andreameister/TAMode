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
#include "sundials_nvector.h"
#include "nvector_serial.h"
#include "reactCode.h"

using namespace std;

static void randomParams (double *pp) {
    for (int ii = 0; ii < 50; ii++) {
        pp[ii] = 10*(double)arc4random() / (double)RAND_MAX;
    }
}




int main() {
    
    double doses[] = {0, 25, 50, 100};
    double pp[50];
    
    TAMout outData;
    outData.pY = (double *) malloc(sizeof(double)*3*NELEMS(doses));
    outData.surf = (double *) malloc(sizeof(double)*3*NELEMS(doses));
    outData.total = (double *) malloc(sizeof(double)*3*NELEMS(doses));
    
    for (size_t trials = 0; trials < 100; trials++) {
        randomParams(pp);
        
        //calcSingleTAMdose (&outData, &params, 10, doses, NELEMS(doses));
        
        double error = calcLew (pp);
        
        cout << error << endl;
        
        
    }
    
    free(outData.pY);
    free(outData.surf);
    free(outData.total);
    
    return 0;
}

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
    for (int ii = 0; ii < 10; ii++) {
        pp[ii] = 10*(double)arc4random() / (double)RAND_MAX;
    }
}




int main() {
    double pp[10];
    
    for (size_t trials = 0; trials < 1000; trials++) {
        randomParams(pp);
        
        double error = calcLew (pp);
        
        if (error < 1.4)
            cout << error << endl;
        
    }
    
    return 0;
}

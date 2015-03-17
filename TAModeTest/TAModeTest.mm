//
//  TAModeTest.m
//  TAModeTest
//
//  Created by Aaron Meyer on 2/17/15.
//  Copyright (c) 2015 Aaron Meyer. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>

#include <stdlib.h>
#include <stdio.h>
#include "cvode.h"
#include <time.h>
#include "TAMode.h"
#include "sundials_nvector.h"
#include "nvector_serial.h"
#include "reactCode.h"

@interface TAModeTest : XCTestCase

@end

@implementation TAModeTest

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
    srand((unsigned int) time(NULL));
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

static struct rates randomParams () {
    double pp[50];
    
    for (int ii = 0; ii < 50; ii++) {
        pp[ii] = 10*(double)arc4random() / (double)RAND_MAX;
    }
    
    struct rates params = Param(pp);
    params.fElse = 0.1;
    params.pD1 = 1;
    
    return params;
}

- (void)testEquilibrium {
    // This test case makes sure that the numbers don't explode at long time in the reaction code.
    N_Vector init = N_VNew_Serial(Nspecies);
    N_Vector dinit = N_VNew_Serial(Nspecies);
    void *cvode_mem;
    
    struct rates params;
    
    for (size_t trials = 0; trials < 100; trials++) {
        cvode_mem = NULL;
        
        params = randomParams();
        
        cvode_mem = initState(init, &params);
        
        AXL_react(0, init, dinit, &params);
        
        XCTAssert(cvode_mem != NULL);
        
        if (cvode_mem != NULL) CVodeFree(&cvode_mem);
        
        XCTAssertEqualWithAccuracy(N_VMaxNorm(dinit), 0, 10);
    }
    
    N_VDestroy(init);
    N_VDestroy(dinit);
}

- (void)testMassConservation {
    // This test case makes sure that the numbers don't explode at long time in the reaction code.
    
    N_Vector init = N_VNew_Serial(Nspecies);
    double totBefore[4];
    void *cvode_mem;
    int flag = 0;
    double t;
    
    struct rates params;
    
    for (size_t trials = 0; trials < 100; trials++) {
        t = 0;
        cvode_mem = NULL;
        
        params = randomParams();
        params.gasCur = params.autocrine;
        params.kDeg = 0;
        params.TAMs[0].expression = 0;
        params.TAMs[1].expression = 0;
        params.TAMs[2].expression = 0;
        
        cvode_mem = initState(init, &params);
        
        XCTAssert(cvode_mem != NULL);
        
        for (int ii = 0; ii < Nspecies; ii++) {
            NV_Ith_S(init,ii) = 10;
        }
        
        totBefore[0] = totCalc(init, &params);
        totBefore[1] = receptorTotCalc(init, &params, AXL);
        totBefore[2] = receptorTotCalc(init, &params, Mer);
        totBefore[3] = receptorTotCalc(init, &params, Tyro);
        
        if (cvode_mem == NULL) continue;
        
        flag = CVodeReInit(cvode_mem, 0.0, init);
        XCTAssertGreaterThanOrEqual(flag, 0);
        
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            continue;
        }
        
        flag = CVode(cvode_mem, autocrineT, init, &t, CV_NORMAL);
        XCTAssertGreaterThanOrEqual(flag, 0);
        if (flag < 0) {
            CVodeFree(&cvode_mem);
            continue;
        }
        
        CVodeFree(&cvode_mem);
        
        XCTAssertEqualWithAccuracy((totBefore[0] - totCalc(init, &params))/totBefore[0], 0, 2E-6);
        XCTAssertEqualWithAccuracy((totBefore[1] - receptorTotCalc(init, &params, AXL))/totBefore[1], 0, 2E-6);
        XCTAssertEqualWithAccuracy((totBefore[2] - receptorTotCalc(init, &params, Mer))/totBefore[2], 0, 2E-6);
        XCTAssertEqualWithAccuracy((totBefore[3] - receptorTotCalc(init, &params, Tyro))/totBefore[3], 0, 2E-6);
    }
    
    N_VDestroy(init);
}

@end

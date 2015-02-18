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
#include "TAMode.h"
#include "sundials_nvector.h"
#include "nvector_serial.h"

@interface TAModeTest : XCTestCase

@end

@implementation TAModeTest

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)testExample {
    // This is an example of a functional test case.
    
    
    double pp[50];
    
    for (int ii = 0; ii < 50; ii++) {
        pp[ii] = 10;
    }
    
    pp[3] = 1;
    pp[1] = 0.1;
    
    void *cvode_mem = NULL;
    struct rates params = Param(pp);
    
    N_Vector init = N_VNew_Serial(Nspecies);
    N_Vector dinit = N_VNew_Serial(Nspecies);
    
    
    
    cvode_mem = initState(init, &params);
    
    int flag = AXL_react(0, init, dinit, &params);
    
    XCTAssert(flag == 0);
    XCTAssert(cvode_mem != NULL);
    
    XCTAssertEqualWithAccuracy(Ith(dinit,4), 0, 0.0001);
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end

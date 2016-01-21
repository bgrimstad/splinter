/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <testingutilities.h>
#include <polynomial.h>
#include "polynomial.h"
#include <polynomialbuilder.h>
#include <testfunctions.h>
#include <testfunction.h>
#include <cmath>
#include "utilities.h"
#include <iostream>
#include "polynomialbuilder.h"

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomialregression][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("Polynomial function" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    double one_eps = 6e-7;
    double two_eps = 6e-7;
    double inf_eps = 6e-7;

    for (auto testFunc : getPolynomialFunctions())
    {
        if (!testFunc->isConstDegree())
        {
            continue;
        }

        auto powers = testFunc->getPowers();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
        CHECK_NOTHROW(compareFunctionValue(testFunc,
                                           [powers](const DataTable &table)
                                           {
                                               Polynomial poly = Polynomial::Builder(table).powers(powers).build();
                                               return new Polynomial(poly);
                                           },
                                           300,  // Number of points to sample at
                                           1337, // Number of points to test against
                                           one_eps, two_eps, inf_eps));
    }
}


TEST_CASE("Polynomial jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]")
{
    double one_eps = 5e-6;
    double two_eps = 5e-6;
    double inf_eps = 5e-5;

    for (auto testFunc : getPolynomialFunctions())
    {
        if (!testFunc->isConstDegree())
        {
            continue;
        }

        auto powers = testFunc->getPowers();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
        CHECK_NOTHROW(compareJacobianValue(testFunc,
                                           [powers](const DataTable &table)
                                           {
                                               Polynomial poly = Polynomial::Builder(table).powers(powers).build();
                                               return new Polynomial(poly);
                                           },
                                           300,  // Number of points to sample at
                                           1337, // Number of points to test against
                                           one_eps, two_eps, inf_eps));
    }
}

//TEST_CASE("Polynomial hessian" COMMON_TEXT, COMMON_TAGS "[hessian]")
//{
//    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
//    {
//        for(auto testFunc : getPolynomialFunctions())
//        {
//            checkHessianSymmetry(testFunc,
//                                 [degree](const DataTable &table)
//                                 {
//                                     return (Function *) new Polynomial(table, degree);
//                                 },
//                                 300,   // Number of points to sample at
//                                 1337); // Number of points to test against
//        }
//    }
//}

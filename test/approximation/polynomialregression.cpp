/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "term.h"
#include <Catch.h>
#include <testingutilities.h>
#include <polynomialregression.h>
#include <testfunctions.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomialregression][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"
#define MAX_DEGREE 4


TEST_CASE("PolynomialRegression function" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
    {
        for(auto testFunc : getPolynomialFunctions())
        {
            double one_eps = 0.1;
            double two_eps = 0.1;
            double inf_eps = 0.1;

            // If the degree of the exact function is less than or equal to the degree
            // of the PolynomialRegression it should fit it exactly.
            if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= degree)
            {
                one_eps = 1e-9;
                two_eps = 1e-9;
                inf_eps = 1e-9;
            }

            CHECK_NOTHROW(compareFunctionValue(testFunc,
                                 [degree](const DataTable &table)
                                 {
                                     return (Approximant *) new PolynomialRegression(table, degree);
                                 },
                                 300,  // Number of points to sample at
                                 1337, // Number of points to test against
                                 one_eps, two_eps, inf_eps));
        }
    }
}

// TODO: Uncomment these when implemented
//TEST_CASE("PolynomialRegression jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]")
//{
//    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
//    {
//        for(auto testFunc : getPolynomialFunctions())
//        {
//            double one_eps = 0.1;
//            double two_eps = 0.1;
//            double inf_eps = 0.1;
//
//            // If the degree of the exact function is less than or equal to the degree
//            // of the PolynomialRegression we are using to approximate it, the approximant
//            // should approximate the function exactly.
//            if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= degree)
//            {
//                one_eps = 1e-5;
//                two_eps = 1e-5;
//                inf_eps = 1e-5;
//            }
//
//            compareJacobianValue(testFunc,
//                                 [degree](const DataTable &table)
//                                 {
//                                     return (Approximant *) new PolynomialRegression(table, degree);
//                                 },
//                                 300,  // Number of points to sample at
//                                 1337, // Number of points to test against
//                                 one_eps, two_eps, inf_eps);
//        }
//    }
//}
//
//TEST_CASE("PolynomialRegression hessian" COMMON_TEXT, COMMON_TAGS "[hessian]")
//{
//    for(int degree = 1; degree <= MAX_DEGREE; ++degree)
//    {
//        for(auto testFunc : getPolynomialFunctions())
//        {
//            checkHessianSymmetry(testFunc,
//                                 [degree](const DataTable &table)
//                                 {
//                                     return (Approximant *) new PolynomialRegression(table, degree);
//                                 },
//                                 300,   // Number of points to sample at
//                                 1337); // Number of points to test against
//        }
//    }
//}

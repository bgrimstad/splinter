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
#include <polynomialregression.h>
#include <testfunctions.h>
#include <testfunction.h>
#include <cmath>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomialregression][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PolynomialRegression function" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    double one_eps = 5e-7;
    double two_eps = 5e-7;
    double inf_eps = 5e-7;

    for (auto testFunc : getPolynomialFunctions())
    {
        if (!testFunc->isConstDegree())
        {
            continue;
        }

        auto degree = testFunc->getConstDegreeInt();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
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


TEST_CASE("PolynomialRegression jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]")
{
    double one_eps = 5e-6;
    double two_eps = 5e-6;
    double inf_eps = 5e-6;

    for (auto testFunc : getPolynomialFunctions())
    {
        if (!testFunc->isConstDegree())
        {
            continue;
        }

        auto degree = testFunc->getConstDegreeInt();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
        CHECK_NOTHROW(compareJacobianValue(testFunc,
                                           [degree](const DataTable &table)
                                           {
                                               return (Approximant *) new PolynomialRegression(table, degree);
                                           },
                                           300,  // Number of points to sample at
                                           1337, // Number of points to test against
                                           one_eps, two_eps, inf_eps));
    }
}

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

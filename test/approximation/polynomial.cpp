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
#include <polynomial2.h>
#include <polynomialbuilder.h>
#include <testfunctions.h>
#include <testfunction.h>
#include <cmath>
#include "utilities.h"
#include <iostream>
#include "polynomial2builder.h"

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

        auto degree = testFunc->getConstDegreeInt();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
        CHECK_NOTHROW(compareFunctionValue(testFunc,
                                           [degree](const DataTable &table)
                                           {
                                               Polynomial poly = Polynomial::Builder(table).degree(degree).build();
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

        auto degree = testFunc->getConstDegreeInt();

        //INFO("Degree: " << degree);
        INFO(testFunc->getFunctionStr());
        CHECK_NOTHROW(compareJacobianValue(testFunc,
                                           [degree](const DataTable &table)
                                           {
                                               Polynomial poly = Polynomial::Builder(table).degree(degree).build();
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


TEST_CASE("Polynomial implementations" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    DataTable dt;

    for (int i = 0; i < 10; i++)
    {
        for (int j = 0; j < 10; j++)
        {
            DenseVector x(2);
            x << i, j;
            double y = 0.1 + 0.2*i - 0.3*i*i + 0.4*i*i*i + 0.5*j - 0.6*i*j - 0.7*j*j;
            dt.addSample(x, y);
        }
    }

    // Original polynomial computes all possible combinations
    std::vector<unsigned int> deg = {3, 2};

    // Coefficient of last monomial with powers x^3*y^3 should be zero
    DenseMatrix pow(8, 2);
    pow << 0, 0,
           1, 0,
           2, 0,
           3, 0,
           0, 1,
           1, 1,
           0, 2,
           3, 3;

    Polynomial poly = Polynomial::Builder(dt).degree(deg).build();
    Polynomial2 poly2 = Polynomial2::Builder(dt).powers(pow).build();
//    Polynomial2 poly2(dt, pow);

    std::cout << "First coefficients" << std::endl;
    std::cout << poly.getCoefficients() << std::endl;

    std::cout << "Sec coefficients" << std::endl;
    std::cout << poly2.getCoefficients() << std::endl;

    std::cout << "Degree: ";
    std::cout << poly2.getDegree() << std::endl;


    for (double xi = 0; xi < 10; xi += 0.01)
    {
        for (double xj = 0; xj < 10; xj += 0.01)
        {
            DenseVector x(2);
            x << xi, xj;
            double y = poly.eval(x);
            double y2 = poly2.eval(x);
            if (!assertNear(y, y2, 1e-10))
                throw Exception("Polynomials not equal!");
        }
    }
}
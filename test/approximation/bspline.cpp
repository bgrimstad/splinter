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

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][bspline][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("Linear BSpline function" COMMON_TEXT, COMMON_TAGS "[bsplinetype::linear][function-value]")
{
    for(auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= 1.0) {
            one_eps = 1e-9;
            two_eps = 1e-9;
            inf_eps = 1e-9;
        }

        compareFunctionValue(testFunc,
                             [](const DataTable &table) {
                                 return (Approximant *) new BSpline(table, BSplineType::LINEAR);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Linear BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::linear][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::LINEAR);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Linear BSpline hessian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::linear][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::LINEAR);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}



TEST_CASE("Quadratic BSpline function" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quadratic][function-value]") {

    for(auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= 2.0) {
            one_eps = 1e-9;
            two_eps = 1e-9;
            inf_eps = 1e-9;
        }

        compareFunctionValue(testFunc,
                             [](const DataTable &table) {
                                 return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quadratic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quadratic][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Quadratic BSpline hessian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quadratic][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}



TEST_CASE("Cubic BSpline function" COMMON_TEXT, COMMON_TAGS "[bsplinetype::cubic][function-value]")
{
    for(auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= 3.0) {
            one_eps = 1e-9;
            two_eps = 1e-9;
            inf_eps = 1e-9;
        }

        compareFunctionValue(testFunc,
                             [](const DataTable &table) {
                                 return (Approximant *) new BSpline(table, BSplineType::CUBIC);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Cubic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::cubic][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::CUBIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Cubic BSpline hessian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::cubic][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::CUBIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}



TEST_CASE("Quartic BSpline function" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quartic][function-value]") {

    for(auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if(testFunc->getF()->isConstDegree() && testFunc->getF()->getConstDegree() <= 4.0) {
            one_eps = 1e-9;
            two_eps = 1e-9;
            inf_eps = 1e-9;
        }

        compareFunctionValue(testFunc,
                             [](const DataTable &table) {
                                 return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quartic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quartic][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Quartic BSpline hessian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quartic][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

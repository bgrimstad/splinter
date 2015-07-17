/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <test_functions.h>
#include <Catch.h>
#include <testingutilities.h>
#include <polynomialregression.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomialregression][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PolynomialRegression function" COMMON_TEXT, COMMON_TAGS "[function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::FunctionValue,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PolynomialRegression jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::Jacobian,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PolynomialRegression hessian" COMMON_TEXT, COMMON_TAGS "[hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::Hessian,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }
}

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
#include <pspline.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][pspline][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PSpline function" COMMON_TEXT, COMMON_TAGS "[function-value][polynomial]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
                      [](const DataTable &table) {
                          return (Approximant *) new PSpline(table);
                      },
                      TestType::FunctionValue,
                      300,  // Number of points to sample at
                      1337, // Number of points to test against
                      one_eps, two_eps, inf_eps);
}

TEST_CASE("PSpline jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian][polynomial]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
                      [](const DataTable &table) {
                          return (Approximant *) new PSpline(table);
                      },
                      TestType::Jacobian,
                      300,  // Number of points to sample at
                      1337, // Number of points to test against
                      one_eps, two_eps, inf_eps);
}

TEST_CASE("PSpline hessian" COMMON_TEXT, COMMON_TAGS "[hessian][polynomial]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
                      [](const DataTable &table) {
                          return (Approximant *) new PSpline(table);
                      },
                      TestType::Hessian,
                      300,  // Number of points to sample at
                      1337, // Number of points to test against
                      one_eps, two_eps, inf_eps);
}

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <test_utils.h>
#include <bspline_builder.h>
#include <test_function.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][pspline]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PSpline function" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.35;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 BSpline bs = BSpline::Builder(table.getDimX(), table.getDimY())
                                         .smoothing(BSpline::Smoothing::PSPLINE)
                                         .alpha(0.03)
                                         .fit(table);
                                 return (Function*) new BSpline(bs);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PSpline function2" COMMON_TEXT, COMMON_TAGS "[function-value2]")
{
    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.8;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 BSpline bs = BSpline::Builder(table.getDimX(), table.getDimY())
                                         .smoothing(BSpline::Smoothing::PSPLINE)
                                         .knotSpacing(BSpline::KnotSpacing::EXPERIMENTAL)
                                         .numBasisFunctions(10)
                                         .alpha(0.01)
                                         .fit(table);
                                 return (Function*) new BSpline(bs);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PSpline jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]")
{
    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 6e-6;
        double two_eps = 6e-6;
        double inf_eps = 6e-5;

        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 BSpline bs = BSpline::Builder(table.getDimX(), table.getDimY())
                                         .smoothing(BSpline::Smoothing::PSPLINE)
                                         .alpha(0.03)
                                         .fit(table);
                                 return (Function*) new BSpline(bs);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

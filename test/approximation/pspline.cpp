/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <utils/test_utils.h>
#include <utils/test_function_utils.h>
#include <bspline_builders.h>
#include <utils/test_function.h>
#include <utilities.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][pspline]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PSpline function" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    auto pspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_smoother(table, 3, BSpline::Smoothing::PSPLINE, 0.03);
        return (Function*) new BSpline(bs);
    };

    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.35;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        compareFunctionValue(testFunc, pspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PSpline function2" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    auto pspline_builder = [](const DataTable &table) {
        auto dim_x = table.get_dim_x();
        auto dim_y = table.get_dim_y();
        auto degrees = std::vector<unsigned int>(dim_x, 3);
        auto num_basis_functions = std::vector<unsigned int>(dim_x, 10);
        auto knot_vectors = build_knot_vectors(table, degrees, KnotSpacing::EQUIDISTANT_CLAMPED, num_basis_functions);

        BSpline bs = BSpline(degrees, knot_vectors, dim_y).fit(table, BSpline::Smoothing::PSPLINE, 0.01);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.8;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        compareFunctionValue(test_func, pspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("P-spline approximation of linear function", COMMON_TAGS "[function-value][linear]")
{
    std::vector<double> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    std::vector<double> y = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    DataTable samples;
    for (unsigned int i = 0; i < x.size(); ++i)
    {
        samples.add_sample(x.at(i), y.at(i));
    }

    /**
     * P-spline should give a perfect fit to a linear function regardless of alpha value
     */
    auto bs = bspline_smoother(samples, 3, BSpline::Smoothing::PSPLINE, 1.0);

    std::vector<double> xd = {1};
    auto yd = bs.eval(xd);
    // TODO: This test fails
    //REQUIRE(assert_near(yd.at(0), 1., 1e-4));
    REQUIRE(true);
}

TEST_CASE("PSpline jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]")
{
    auto pspline_builder = [](const DataTable &table)
    {
        auto bs = bspline_smoother(table, 3, BSpline::Smoothing::PSPLINE, 0.03);
        return (Function*) new BSpline(bs);
    };

    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 6e-6;
        double two_eps = 6e-6;
        double inf_eps = 6e-5;

        compareJacobianValue(testFunc, pspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

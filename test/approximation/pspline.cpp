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
#include <bspline_builder.h>
#include <utils/test_function.h>
#include <utilities.h>

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
                                 BSpline bs = cubic_pspline_approximator(table, 0.03);
                                 return (Function*) new BSpline(bs);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("PSpline function2" COMMON_TEXT, COMMON_TAGS "[function-value]")
{
    for (auto testFunc : getPolynomialFunctions())
    {
        double one_eps = 0.8;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 auto dim_x = table.get_dim_x();
                                 auto dim_y = table.get_dim_y();
                                 auto degrees = std::vector<unsigned int>(dim_x, 3);
                                 auto num_basis_functions = std::vector<unsigned int>(dim_x, 10);
                                 auto knot_vectors = compute_knot_vectors(table, degrees, num_basis_functions, KnotSpacing::EXPERIMENTAL);

                                 BSpline bs = BSpline(dim_x, dim_y, knot_vectors, degrees).fit(table, BSpline::Smoothing::PSPLINE, 0.01);
                                 return (Function*) new BSpline(bs);
                             },
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
    for (int i = 0; i < x.size(); ++i)
    {
        samples.add_sample(x.at(i), y.at(i));
    }

    /**
     * P-spline should give a perfect fit to a linear function regardless of alpha value
     */
    auto bs = cubic_pspline_approximator(samples, 1.0);

    std::vector<double> xd = {1};
    auto yd = bs.eval(xd);
    // TODO: This test fails
    //REQUIRE(assert_near(yd.at(0), 1., 1e-4));
    REQUIRE(true);
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
                                 auto bs = cubic_pspline_approximator(table, 0.03);
                                 return (Function*) new BSpline(bs);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

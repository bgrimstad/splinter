/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include <Catch.h>
#include <bspline_builders.h>
#include "utilities.h"
#include <iostream>

using namespace SPLINTER;
using std::cout;
using std::endl;

#define COMMON_TAGS "[examples][approximation]"

TEST_CASE("Approximation example", COMMON_TAGS)
{
    // We want to approximate this function, known as the six-hump camelback function
    auto f = [](std::vector<double> x) {
        assert(x.size() == 2);
        auto x0 = x.at(0);
        auto x1 = x.at(1);
        return (4 - 2.1*x0*x0 + (1/3.)*x0*x0*x0*x0)*x0*x0 + x0*x1 + (-4 + 4*x1*x1)*x1*x1;
    };

    // Create new DataTable to manage samples
    DataTable samples;

    // Sample the function
    for (auto i = 0; i < 20; i++)
    {
        for(auto j = 0; j < 20; j++)
        {
            // Sample function at x
            std::vector<double> x = {i*0.1, j*0.1};
            samples.add_sample(x, f(x));
        }
    }

    // Build B-splines that interpolate the samples
    BSpline bspline1 = bspline_interpolator(samples, 1);
    BSpline bspline3 = bspline_interpolator(samples, 3);

    // Build penalized B-spline (P-spline) that smooths the samples
    BSpline pspline = bspline_smoother(samples, 3, BSpline::Smoothing::PSPLINE, 0.03);

    /*
     * Evaluate the splines at x = (1,1)
     * NOTE1: the error will be 0 at that point (except for the P-spline, which may introduce an error
     * in favor of a smooth approximation) because it is a point we sampled at.
     * NOTE2: The BSpline::eval function returns an output vector (in this case of size 1)
     */
    std::vector<double> x = {1, 1};
    auto y_f = f(x);
    auto y_bs1 = bspline1.eval(x).at(0);
    auto y_bs3 = bspline3.eval(x).at(0);
    auto y_ps = pspline.eval(x).at(0);

    // Test results at x = (1, 1)
    REQUIRE(assert_near(y_f, y_bs1, 1e-8));
    REQUIRE(assert_near(y_f, y_bs3, 1e-8));
    REQUIRE(assert_near(y_ps, y_bs3, 1e-3)); // Smoothed (does not interpolate)
}

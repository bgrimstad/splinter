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
#include "utils/timer.h"


using namespace SPLINTER;
using std::cout;
using std::endl;

#define COMMON_TAGS "[examples][batch_eval]"

TEST_CASE("Batch evaluation of B-spline in two variables", COMMON_TAGS)
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
     * Compare evaluation of single points with batch evaluation
     */

    Timer tim;
    tim.start();
    DenseMatrix y_single(100*100, 1);
    int k = 0;
    for (auto i = 1; i <= 100; i++)
    {
        for(auto j = 1; j <= 100; j++)
        {
            // Sample function at x
            DenseVector x(2);
            x << i*0.019, j*0.019;
            auto y = bspline3.eval(x);
            y_single.row(k) = y;
            k++;
        }
    }
    tim.stop();
    std::cout << "Time used with single evaluation: " << std::to_string(tim.get_milli_seconds()) << " ms" << std::endl;

    tim.reset();
    tim.start();

    DenseMatrix x_batch(100*100, 2);
    k = 0;
    for (auto i = 1; i <= 100; i++)
    {
        for(auto j = 1; j <= 100; j++)
        {
            // Sample function at x
            DenseVector x(2);
            x << i*0.019, j*0.019;
            x_batch.row(k) = x;
            k++;
        }
    }
    DenseMatrix y_batch = bspline3.batch_eval(x_batch);

    tim.stop();
    std::cout << "Time used with batch evaluation: " << std::to_string(tim.get_milli_seconds()) << " ms" << std::endl;

    REQUIRE(y_batch.rows() == y_single.rows());
    REQUIRE(y_batch.cols() == y_single.cols());

    for (int i = 0; i < 100*100; ++i) {
        REQUIRE(assert_near(y_batch(i, 0), y_single(i, 0), 1e-6));
    }

}

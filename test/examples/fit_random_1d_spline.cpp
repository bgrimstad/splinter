/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * This example was provided by IuriiSorokin (https://github.com/IuriiSorokin)
 */

#include <bspline.h>
#include <data_table.h>
#include <utilities.h>
#include <Catch.h>
#include <cmath>
#include <random>

using namespace SPLINTER;

#define COMMON_TAGS "[examples][approximation]"

TEST_CASE("Random 1-D spline example", COMMON_TAGS)
{
    auto random_engine = std::default_random_engine();

    // Function, the data points will follow
    const auto func = []( double x ){ return std::sin(x); };

    // Range
    const double x_min = 0;
    const double x_max = 2*M_PI;

    // Scattering of the points w.r.t. the function
    const double sigma = 1;
    auto scattering_distr = std::normal_distribution<double>(0, sigma);

    // Generate data
    const size_t n_points = 150;
    DataTable data_table;
    auto x_distr = std::uniform_real_distribution<double>(x_min, x_max);
    for( size_t i_point = 0; i_point < n_points; ++i_point )
    {
        const double x = x_distr( random_engine );
        const double y = func(x) + scattering_distr(random_engine);
        data_table.add_sample(x, y);
    }

    // Define and fit spline
    const unsigned int y_dim = 1;
    const std::vector<unsigned int> degrees = {3};
    const std::vector<std::vector<double>> knots = { {x_min, x_min, x_min, x_min, x_max, x_max, x_max, x_max } };
    BSpline spline(degrees, knots, y_dim);
    spline.fit(data_table);

    auto y1 = spline.eval(std::vector<double>({0})).at(0);
    auto y2 = spline.eval(std::vector<double>({0.5})).at(0);
    auto y3 = spline.eval(std::vector<double>({M_PI})).at(0);
    auto y4 = spline.eval(std::vector<double>({2*M_PI})).at(0);

    CHECK(assert_near(y1, 0.08363589561));
    CHECK(assert_near(y2, 0.6305184164));
    CHECK(assert_near(y3, 0.0237115102));
    CHECK(assert_near(y4, -0.08151616982));
}

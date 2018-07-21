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
#include <utilities.h>
#include <bspline_builders.h>

using namespace SPLINTER;

#define COMMON_TAGS "[approximation][bspline][wls]"


// TESTING LINEAR B-SPLINES
TEST_CASE("Weighted least squares with linear B-spline", COMMON_TAGS "[bsplinetype::linear]")
{
    DataTable data;

    unsigned int n = 10;
    auto x = linspace(1, 5, n);
    auto y = linspace(1, 5, n);

    // Set outlying sample
    unsigned int m = 2;
    auto y_true = y.at(m);
    y.at(m) = 1e5;

    for (unsigned int i = 0; i < n; ++i)
        data.add_sample(x.at(i), y.at(i));

    // Set weights - zero weight on outlying sample
    auto weights = std::vector<double>(n, 1);
    weights.at(m) = 0.0;

    // Fit B-spline using WLS
    // NOTE: It does not make sense to use WLS when interpolating all data points - in this case the linear system may
    // become ill-conditioned if some of the weights are set to zero.
    auto degrees = std::vector<unsigned int>(1, 1);
    auto num_basis_functions = std::vector<unsigned int>(1, int(n/2));
    auto knot_vectors = build_knot_vectors(data, degrees, KnotSpacing::EQUIDISTANT_CLAMPED, num_basis_functions);

    auto bs = BSpline(degrees, knot_vectors).fit(data, BSpline::Smoothing::NONE, 0., weights);
    std::vector<double> x_eval = {x.at(m)};
    auto y_eval = bs.eval(x_eval).at(0);
    REQUIRE(assert_near(y_eval, y_true));
}

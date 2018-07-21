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
#include <utilities.h>

using namespace SPLINTER;

#define COMMON_TAGS "[general][bsplinebuilder]"
#define COMMON_TEXT "BSplineBuilder "

TEST_CASE(COMMON_TEXT "multivariate output", COMMON_TAGS "[construction]")
{
    auto table_y1 = DataTable();
    auto table_y2 = DataTable();
    auto table_ys = DataTable();

    auto x_vec = linspace(0, 10, 11);

    auto f1 = [](double x) { return 2*x*x - x; };
    auto f2 = [](double x) { return x*x*x; };

    for (auto xi : x_vec)
    {
        table_y1.add_sample(xi, f1(xi));
        table_y2.add_sample(xi, f2(xi));
        table_ys.add_sample(xi, {f1(xi), f2(xi)});
    }

    unsigned int degree = 3;
    auto bs_y1 = bspline_interpolator(table_y1, degree);
    auto bs_y2 = bspline_interpolator(table_y2, degree);
    auto bs_ys = bspline_interpolator(table_ys, degree);

    auto x_vec_test = linspace(-10, 20, 100);

    for (auto xi : x_vec_test)
    {
        auto xi_ = std::vector<double>({xi});
        auto y1 = bs_y1.eval(xi_);
        auto y2 = bs_y2.eval(xi_);
        auto ys = bs_ys.eval(xi_);

        assert(assert_near(y1.at(0), ys.at(0)));
        assert(assert_near(y2.at(0), ys.at(1)));
    }
}

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <bspline_builder.h>
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
        table_y1.addSample(xi, f1(xi));
        table_y2.addSample(xi, f2(xi));
        table_ys.addSample(xi, {f1(xi), f2(xi)});
    }

    auto bs_y1 = BSpline::Builder(table_y1).degree(3).build();
    auto bs_y2 = BSpline::Builder(table_y2).degree(3).build();
    auto bs_ys = BSpline::Builder(table_ys).degree(3).build();

    auto x_vec_test = linspace(-10, 20, 100);

    for (auto xi : x_vec_test)
    {
        auto xi_ = std::vector<double>({xi});
        auto y1 = bs_y1.eval(xi_);
        auto y2 = bs_y2.eval(xi_);
        auto ys = bs_ys.eval(xi_);

        assert(assertNear(y1.at(0), ys.at(0)));
        assert(assertNear(y2.at(0), ys.at(1)));
    }
}

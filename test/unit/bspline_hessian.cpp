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

#define COMMON_TAGS "[unit][bspline][hessian]"
#define COMMON_TEXT " unit test"

TEST_CASE("hessian" COMMON_TEXT, COMMON_TAGS)
{
    /*
     * Test functions
     */
    auto f1 = [](double x1, double x2) { return 1 + 2*x1 + 3*x2 + 4*x1*x2 + 5*x1*x1 + 6*x2*x2; };
    auto f2 = [](double x1, double x2) { return 6 + 5*x1 + 4*x2 + 3*x1*x2 + 2*x1*x1 + 1*x2*x2; };

    /*
     * Expected Hessian at (1, 1)
     * |10  4| and |4 3|
     * |4  12|     |3 2|
     */
    DenseMatrix f1_hessian(2, 2);
    f1_hessian << 10, 4, 4, 12;
    DenseMatrix f2_hessian(2, 2);
    f2_hessian << 4, 3, 3, 2;
    std::vector<DenseMatrix> hessian_true = {f1_hessian, f2_hessian};

    /*
     * Sample and create B-spline
     */
    auto x1_vec = linspace(-1, 1, 10);
    auto x2_vec = linspace(0, 2, 10);

    DataTable samples;
    for (auto x1_i : x1_vec)
        for (auto x2_i : x2_vec)
            samples.add_sample(std::vector<double>({x1_i, x2_i}),
                               std::vector<double>({f1(x1_i, x2_i), f2(x1_i, x2_i)}));

    BSpline bs = bspline_interpolator(samples, 2);

    /*
     * Evaluate Hessian and compare to true Hessian
     */
    auto hessian = bs.eval_hessian({1, 1});

    for (size_t i = 0; i < hessian.size(); ++i) {
        for (size_t j = 0; j < hessian.at(i).size(); ++j) {
            for (size_t k = 0; k < hessian.at(i).at(j).size(); ++k) {
                double Hijk_calc = hessian.at(i).at(j).at(k);
                double Hijk_true = hessian_true.at(i)(j, k);
                REQUIRE(assert_near(Hijk_calc, Hijk_true));
            }
        }
    }
}

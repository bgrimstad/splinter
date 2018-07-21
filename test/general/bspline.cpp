/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <utils/bspline_test_utils.h>
#include <utilities.h>

using namespace SPLINTER;

#define COMMON_TAGS "[general][bspline]"
#define COMMON_TEXT "BSpline "


TEST_CASE(COMMON_TEXT "construction multivariate", COMMON_TAGS "[construction]")
{
    // Build a tensor product B-spline f : R -> R
    unsigned int num_knots = 11;
    std::vector<std::vector<double>> knots = {linspace(0, 10, num_knots)};
    std::vector<unsigned int> deg = {3};

    unsigned int num_cp = num_knots - deg.at(0) - 1;
    auto cp1_vec = linspace(0, 6, num_cp);
    std::vector<std::vector<double>> cp1;
    for (auto cpi : cp1_vec)
        cp1.push_back({cpi});

    auto cp2_vec = linspace(6, 13, num_cp);
    std::vector<std::vector<double>> cp2;
    for (auto cpi : cp2_vec)
        cp2.push_back({cpi});

    // Build two 1-D B-splines
    auto bs1 = BSpline(deg, knots, cp1);
    auto bs2 = BSpline(deg, knots, cp2);

    // Build a B-spline f : R -> R^2
    std::vector<std::vector<double>> cp3;
    for (unsigned int i = 0; i < num_cp; ++i)
        cp3.push_back({cp1_vec.at(i), cp2_vec.at(i)});

    auto bs3 = BSpline(deg, knots, cp3);
}

TEST_CASE(COMMON_TEXT "construction throws", COMMON_TAGS "[construction]")
{
    // Build a tensor product B-spline with three input variables
    std::vector<std::vector<double>> kv = {linspace(-10, 10, 21),
                                           linspace(-10, 10, 21),
                                           linspace(-10, 10, 21)};

    std::vector<unsigned int> deg = {1, 2, 3};

    // Building a B-spline without specifying control points (just to get the number of basis functions)
    BSpline bs = BSpline(deg, kv);

    // Create vector of 1-D control points
    // For m outputs and l control points, the constructor expects a vector containing l vectors of size m
    // Here, we construct a vector containing m vectors of size l,
    // which should cause the constructor to throw an exception
    auto cp_vec = std::vector<std::vector<double>>(1, linspace(0, 100, bs.get_num_basis_functions()));

    // Expecting constructor to throw
    REQUIRE_THROWS(BSpline(deg, kv, cp_vec));
}

TEST_CASE(COMMON_TEXT "domain reduction", COMMON_TAGS "[subdivision]")
{
    REQUIRE(domain_reduction_test1());
}

TEST_CASE(COMMON_TEXT "recursive subdivision", COMMON_TAGS "[subdivision]")
{
    // TODO: The current code for comparing BSplines require identical bounds which fails in this test.
    //REQUIRE(run_recursive_domain_reduction_test());
}

TEST_CASE(COMMON_TEXT "knot insertion", COMMON_TAGS "[knotinsertion]")
{
    REQUIRE(test_knot_insertion());
}

TEST_CASE(COMMON_TEXT "knot averages", COMMON_TAGS "[knotaverages]")
{
    // Build a tensor product B-spline with three input variables
    std::vector<std::vector<double>> kv1 = {linspace(-10, 10, 21),
                                            linspace(-10, 10, 21),
                                            linspace(-10, 10, 21)};

    std::vector<unsigned int> deg1 = {1, 2, 3};

    // Building a B-spline without specifying control points (just to get the number of basis functions)
    BSpline bs1 = BSpline(deg1, kv1);

    // Create vector of 1-D control points
    auto cp1_vec = linspace(0, 100, bs1.get_num_basis_functions());
    std::vector<std::vector<double>> cp1;
    for (auto cpi : cp1_vec)
        cp1.push_back({cpi});

    // Building B-spline
    bs1 = BSpline(deg1, kv1, cp1);

    // Get control points and knot averages
    auto cp1_mat = bs1.get_control_points();
    auto mu1_mat = bs1.get_knot_averages();

    // Build a B-spline with multidimensional control points P_i = (mu_i, cp_i), where {mu_i} are the knot averages
    DenseMatrix cp1_new = DenseMatrix::Zero(cp1_mat.rows(), mu1_mat.cols() + 1);
    cp1_new.block(0, 0, mu1_mat.rows(), mu1_mat.cols()) = mu1_mat;
    cp1_new.block(0, mu1_mat.cols(), cp1_mat.rows(), cp1_mat.cols()) = cp1_mat;

    BSpline bs1_new = BSpline(deg1, kv1, eig_to_std_mat(cp1_new));

    /*
     * The new B-spline should evaluate to (x, f(x)) for any x in the B-spline support
     * E.g. for x = (0, 0, 0) we expect y = (0, 0, 0, 50), where 50 is the old B-spline value at x = (0, 0, 0)
     */
    auto x = std::vector<double>({0, 0, 0});
    auto y = bs1_new.eval(x);

    for (size_t i = 0; i < x.size(); ++i)
        REQUIRE(assert_near(x.at(i), y.at(i)));

    REQUIRE(assert_near(50., y.back()));

    /*
     * Testing for many points in the range of the knot averages
     * NOTE: The range of the knot averages is [t_p, t_n] for a knot vector [t_0, ..., t_{n+p}].
     * In this example, the range of the knot averages are [-9, 9], [-8, 8], and [-7, 7].
     */
    bool test = true;
    for (auto x0 : linspace(-9, 9, 50))
    for (auto x1 : linspace(-8, 8, 50))
    for (auto x2 : linspace(-7, 7, 50))
    {
        auto x = std::vector<double>({x0, x1, x2});
        auto y = bs1_new.eval(x);

        if (!assert_near(x0, y.at(0)))
        {
            test = false;
            break;
        }

        if (!assert_near(x1, y.at(1)))
        {
            test = false;
            break;
        }

        if (!assert_near(x2, y.at(2)))
        {
            test = false;
            break;
        }
    }

    REQUIRE(test);
}
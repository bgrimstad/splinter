/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <bspline_basis_1d.h>

using namespace SPLINTER;

#define COMMON_TAGS "[unit][bsplinebasis1d]"
#define COMMON_TEXT " unit test"

TEST_CASE("support_hack" COMMON_TEXT, COMMON_TAGS)
{
    std::vector<double> knots = {1, 1, 1, 2.1, 3.1, 4, 4, 4};
    BSplineBasis1D bb(2, knots);

    double x = 4;
    x = bb.support_hack(x);
    REQUIRE(x < 4);
}

TEST_CASE("index_supported_basis_functions" COMMON_TEXT, COMMON_TAGS)
{
    std::vector<double> knots1 = {1, 1, 1, 2.1, 2.5, 2.5, 2.8, 2.8, 2.8, 3.1, 4, 4, 4};
    BSplineBasis1D bb1(2, knots1);

    {
        double x = 1;
        std::vector<int> ref = {0, 1, 2};
        auto sup = bb1.index_supported_basis_functions(x);

        for (unsigned int i = 0; i < ref.size(); ++i)
            REQUIRE(ref.at(i) == sup.at(i));
    }

    {
        double x = 2.5;
        std::vector<int> ref = {3, 4, 5};
        auto sup = bb1.index_supported_basis_functions(x);

        for (unsigned int i = 0; i < ref.size(); ++i)
            REQUIRE(ref.at(i) == sup.at(i));
    }

    {
        double x = 3.999;
        std::vector<int> ref = {7, 8, 9};
        auto sup = bb1.index_supported_basis_functions(x);

        for (unsigned int i = 0; i < ref.size(); ++i)
            REQUIRE(ref.at(i) == sup.at(i));
    }

    std::vector<double> knots2 = {0, 1, 2.1, 2.5, 2.8, 2.8, 3.1, 4, 4};
    BSplineBasis1D bb2(1, knots2);

    {
        double x = 0;
        std::vector<int> ref = {0};
        auto sup = bb2.index_supported_basis_functions(x);

        for (unsigned int i = 0; i < ref.size(); ++i)
            REQUIRE(ref.at(i) == sup.at(i));

    }

    {
        double x = 2.499999999;
        std::vector<int> ref = {1, 2};
        auto sup = bb2.index_supported_basis_functions(x);

        for (unsigned int i = 0; i < ref.size(); ++i)
            REQUIRE(ref.at(i) == sup.at(i));

    }
}
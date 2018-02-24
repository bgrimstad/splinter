/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <knot_builders.h>
#include <knot_vector.h>
#include <iostream>

using namespace SPLINTER;

#define COMMON_TAGS "[unit][knots]"
#define COMMON_TEXT " unit test"


TEST_CASE("is_regular" COMMON_TEXT, COMMON_TAGS)
{
    auto knots0 = KnotVector({1, 2.1, 3.1, 4});
    auto knots1 = KnotVector({1, 1, 2.1, 3.1, 4, 4});
    auto knots2 = KnotVector({1, 1, 1, 2.1, 3.1, 4, 4, 4});
    auto knots3 = KnotVector({1, 1, 1, 1, 2.1, 3.1, 4, 4, 4, 4});
    auto knots4 = KnotVector({1, 1, 2.0, 2.1, 4, 4});

    REQUIRE(knots0.is_regular(0));
    REQUIRE(knots1.is_regular(1));
    REQUIRE(!knots2.is_regular(1));
    REQUIRE(knots2.is_regular(2));
    REQUIRE(!knots3.is_regular(2));
    REQUIRE(knots3.is_regular(3));
    REQUIRE(knots3.is_regular(4));
    REQUIRE(knots4.is_regular(1));
}

TEST_CASE("is_clamped" COMMON_TEXT, COMMON_TAGS)
{
    auto knots1 = KnotVector({1, 1, 2.1, 3.1, 4, 4});
    auto knots2 = KnotVector({1, 1, 1, 2.1, 3.1, 4, 4, 4});
    auto knots3 = KnotVector({1, 1, 1, 1, 2.1, 3.1, 4, 4, 4, 4});

    REQUIRE(knots1.is_clamped(1));
    REQUIRE(knots2.is_clamped(2));
    REQUIRE(knots3.is_clamped(3));
}

TEST_CASE("is_refinement" COMMON_TEXT, COMMON_TAGS)
{
    auto knots1 = KnotVector({1, 1, 1, 2.1, 3.1, 4, 4, 4});
    auto knots2 = KnotVector({1, 1, 1, 2.1, 2.5, 3.1, 4, 4, 4});

    REQUIRE(knots1.is_refinement(knots2));
}

TEST_CASE("knot_vector_expanded_equidistant" COMMON_TEXT, COMMON_TAGS)
{
    std::vector<double> values = {1, 1, 1, 2.1, 3.1, 4, 4, 4, 3.1, 2.1, 2.2};
    unsigned int degree = 3;
    unsigned int num_basis_functions = 1;
    auto knots = knot_vector_expanded_equidistant(values, degree, num_basis_functions);
    std::vector<double> correct_knots = {0.7, 1.6, 2.5, 3.4, 4.3};

    REQUIRE(knots.size() == correct_knots.size());
    bool flag = true;
    for (unsigned int i = 0; i < knots.size(); ++i)
        if (std::fabs(knots.at(i) - correct_knots.at(i)) > 1e-6)
            flag = false;
    REQUIRE(flag);
}

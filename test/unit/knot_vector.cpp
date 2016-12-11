/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <knot_vector.h>
#include <iostream>
#include <algorithm>

using namespace SPLINTER;

#define COMMON_TAGS "[unit][knots]"
#define COMMON_TEXT " unit test"


TEST_CASE("knot_vector_initialization" COMMON_TEXT, COMMON_TAGS)
{
    std::vector<double> knots0 = {1, 2.1, 3.1, 4};
    std::vector<double> knots1 = {1, 2.1, 3.1, 4, 5};
    std::vector<double> knots2 = {1, 1, 2.1, 3.1, 4, 4};
    std::vector<double> knots3 = {1, 1, 1, 2.1, 3.1, 4, 4, 4};
    std::vector<double> knots4 = {1, 1, 1, 1, 2.1, 3.1, 4, 4, 4, 4};

    auto knot_vector0 = KnotVector(knots0);
    auto knot_vector1 = KnotVector(knots1);
    auto knot_vector2 = KnotVector(knots2);
    auto knot_vector3 = KnotVector(knots3);
    auto knot_vector4 = KnotVector(knots4);

    REQUIRE(knot_vector0 == knot_vector0);
    REQUIRE(knot_vector0 != knot_vector1);
}

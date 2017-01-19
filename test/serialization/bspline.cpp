/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <data_table.h>
#include <bspline_builder.h>
#include "test_utils.h"

using namespace SPLINTER;


#define COMMON_TAGS "[serialization][bspline]"


TEST_CASE("BSpline can be saved and loaded", COMMON_TAGS)
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = multi_linspace(dim, std::pow(300, 1.0 / dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.bspline";

    SECTION("Linear BSpline")
    {
        BSpline bspline = BSpline::Builder(table).degree(1).fit();
        bspline.save(fileName);
        BSpline loadedBSpline(fileName);
        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Quadratic BSpline")
    {
        BSpline bspline = BSpline::Builder(table).degree(2).fit();
        bspline.save(fileName);
        BSpline loadedBSpline(fileName);
        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Cubic BSpline")
    {
        BSpline bspline = BSpline::Builder(table).degree(3).fit();
        bspline.save(fileName);
        BSpline loadedBSpline(fileName);
        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Quartic BSpline")
    {
        BSpline bspline = BSpline::Builder(table).degree(4).fit();
        bspline.save(fileName);
        BSpline loadedBSpline(fileName);
        REQUIRE(bspline == loadedBSpline);
    }

    remove(fileName);
}

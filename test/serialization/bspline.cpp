/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <datatable.h>
#include <bspline.h>
#include "testingutilities.h"

#include "term.h"

using namespace SPLINTER;


TEST_CASE("BSplines can be saved and loaded", "[serialize][bspline]")
{
    auto func = getTestFunction(1, 1);
    auto dim = func->getNumVariables();
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(100, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.bspline";

    SECTION("Linear BSpline") {
        BSpline bspline(table, BSplineType::LINEAR);

        bspline.save(fileName);
        BSpline loadedBSpline(fileName);

        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Quadratic BSpline") {
        BSpline bspline(table, BSplineType::QUADRATIC);

        bspline.save(fileName);
        BSpline loadedBSpline(fileName);

        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Cubic BSpline") {
        BSpline bspline(table, BSplineType::CUBIC);

        bspline.save(fileName);
        BSpline loadedBSpline(fileName);

        REQUIRE(bspline == loadedBSpline);
    }

    SECTION("Quartic BSpline") {
        BSpline bspline(table, BSplineType::QUARTIC);

        bspline.save(fileName);
        BSpline loadedBSpline(fileName);

        REQUIRE(bspline == loadedBSpline);
    }

    remove(fileName);
}

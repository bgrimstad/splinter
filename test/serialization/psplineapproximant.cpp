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
#include <psplineapproximant.h>
#include "testingutilities.h"

using namespace SPLINTER;


TEST_CASE("PSplines can be saved and loaded", "[serialization][pspline]")
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(300, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.pspline";

    SECTION("PSpline with default lambda") {
        PSplineApproximant pspline(table);

        pspline.save(fileName);
        PSplineApproximant loadedPSpline(fileName);

        REQUIRE(pspline == loadedPSpline);
    }

    SECTION("PSpline with non-default lambda") {
        PSplineApproximant pspline(table, 0.02);

        pspline.save(fileName);
        PSplineApproximant loadedPSpline(fileName);

        REQUIRE(pspline == loadedPSpline);
    }

    remove(fileName);
}

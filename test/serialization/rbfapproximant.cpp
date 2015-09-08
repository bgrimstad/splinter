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
#include <rbfapproximant.h>
#include "testingutilities.h"

using namespace SPLINTER;


TEST_CASE("RadialBasisFunction can be saved and loaded", "[serialization][rbf][radialbasisfunction]")
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(300, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.rbf";

    SECTION("Gaussian") {
        RBFApproximant rbf(table, RBFType::GAUSSIAN);

        rbf.save(fileName);
        RBFApproximant loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Inverse multiquadric") {
        RBFApproximant rbf(table, RBFType::INVERSE_MULTIQUADRIC);

        rbf.save(fileName);
        RBFApproximant loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Inverse quadric") {
        RBFApproximant rbf(table, RBFType::INVERSE_QUADRIC);

        rbf.save(fileName);
        RBFApproximant loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Multiquadric") {
        RBFApproximant rbf(table, RBFType::MULTIQUADRIC);

        rbf.save(fileName);
        RBFApproximant loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Thin plate spline") {
        RBFApproximant rbf(table, RBFType::THIN_PLATE_SPLINE);

        rbf.save(fileName);
        RBFApproximant loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    remove(fileName);
}

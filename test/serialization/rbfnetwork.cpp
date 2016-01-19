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
#include "rbfnetwork.h"
#include "rbfbuilder.h"
#include "testingutilities.h"

using namespace SPLINTER;


#define COMMON_TAGS "[serialization][rbf][radialbasisfunction]"


TEST_CASE("RadialBasisFunction can be saved and loaded", COMMON_TAGS)
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(300, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.rbf";

    SECTION("Gaussian") {
        RBFNetwork rbf = RBFNetwork::Builder(table).type(RBFType::GAUSSIAN).build();

        rbf.save(fileName);
        RBFNetwork loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Inverse multiquadric") {
        RBFNetwork rbf = RBFNetwork::Builder(table).type(RBFType::INVERSE_MULTIQUADRIC).build();

        rbf.save(fileName);
        RBFNetwork loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Inverse quadric") {
        RBFNetwork rbf = RBFNetwork::Builder(table).type(RBFType::INVERSE_QUADRIC).build();

        rbf.save(fileName);
        RBFNetwork loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Multiquadric") {
        RBFNetwork rbf = RBFNetwork::Builder(table).type(RBFType::MULTIQUADRIC).build();

        rbf.save(fileName);
        RBFNetwork loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    SECTION("Thin plate spline") {
        RBFNetwork rbf = RBFNetwork::Builder(table).type(RBFType::THIN_PLATE_SPLINE).build();

        rbf.save(fileName);
        RBFNetwork loadedRBF(fileName);

        REQUIRE(rbf == loadedRBF);
    }

    remove(fileName);
}

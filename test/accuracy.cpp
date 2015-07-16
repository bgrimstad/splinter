/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <test_functions.h>
#include <Catch.h>
#include <testingutilities.h>

using namespace SPLINTER;

TEST_CASE("Quartic BSpline accuracy test with nice functions", "[bspline][nice]") {
    auto niceFunctions = getNiceTestFunctions();
    for(auto &exact : niceFunctions) {
        INFO("Function: " << *(exact->getF()));

        auto dim = exact->getNumVariables();
        REQUIRE(dim > 0);
        auto samplePoints = linspace(dim, -5, 5, std::pow(500, 1.0/dim));
        auto evalPoints = linspace(dim, -4.95, 4.95, std::pow(1337, 1.0/dim));

        DataTable table = sample(exact, samplePoints);

        BSpline b(table, BSplineType::CUBIC);

        compareFunctions(*exact, b, evalPoints);
    }
}

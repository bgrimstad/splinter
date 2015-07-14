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

TEST_CASE("BSpline accuracy test", "[bspline]") {
//    Var x(0);
//    Var y(1);
//    Var z(2);
//    auto f = z * ((x^2) + (y^2));
//
//    auto dim = 3;
//
//    auto fewSamplePoints = linspace(dim, -5, 5, 7);
//    auto manySamplePoints = linspace(dim, -5, 5, 10);
//    auto evalPoints = linspace(dim, -4.95, 4.95, 11);
//
//    TestFunction exact(dim, f);
//    DataTable fewSamples = sample(exact, fewSamplePoints);
//    DataTable manySamples = sample(exact, manySamplePoints);

    auto exact = testFunctions.at(0);
    auto dim = exact->getNumVariables();
    auto samplePoints = linspace(dim, -5, 5, 10);
    auto evalPoints = linspace(dim, -4.95, 4.95, 11);

    DataTable table = sample(*exact, samplePoints);

    BSpline b(table, BSplineType::CUBIC);

    INFO("Function: " << *(exact->getF()));

    compareFunctions(*exact, b, evalPoints);
//    for(auto &exact : testFunctions) {
//        auto dim = exact->getNumVariables();
//        auto samplePoints = linspace(dim, -5, 5, 10);
//        auto evalPoints = linspace(dim, -4.95, 4.95, 11);

//        DataTable table = sample(*exact, samplePoints);

//        BSpline b(table, BSplineType::CUBIC);

//        INFO("Function: " << *(exact->getF()));

//        compareFunctions(*exact, b, evalPoints);
//    }
}

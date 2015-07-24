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
#include "polynomialregression.h"
#include "testingutilities.h"

using namespace SPLINTER;


TEST_CASE("PolynomialRegression can be saved and loaded", "[serialization][polynomialregression]")
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(300, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.polyfit";

    SECTION("PolynomialRegression of degree 1") {
        PolynomialRegression polyfit(table, 1);

        polyfit.save(fileName);
        PolynomialRegression loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    SECTION("PolynomialRegression of degree 5") {
        PolynomialRegression polyfit(table, 4);

        polyfit.save(fileName);
        PolynomialRegression loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    if(dim > 1) {
        SECTION("PolynomialRegression of differing degrees in each dimension") {
            auto degrees = std::vector<unsigned int>(dim);
            for(unsigned int i = 0; i < dim; i++) {
                degrees.at(i) = i + 1;
            }
            PolynomialRegression polyfit(table, degrees);

            polyfit.save(fileName);
            PolynomialRegression loadedPolyfit(fileName);

            REQUIRE(polyfit == loadedPolyfit);
        }
    }

    remove(fileName);
}

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
#include "polynomial.h"
#include "polynomialbuilder.h"
#include "testingutilities.h"

using namespace SPLINTER;


#define COMMON_TAGS "[serialization][polynomial]"


TEST_CASE("Polynomial can be saved and loaded", COMMON_TAGS)
{
    unsigned int dim = 2;
    auto func = getTestFunction(dim, 1);
    // Don't sample too fine, this test isn't supposed to test the speed
    auto points = linspace(dim, std::pow(300, 1.0/dim));
    DataTable table = sample(func, points);

    const char *fileName = "test.polyfit";

    SECTION("Polynomial of degree 1")
    {
        Polynomial polyfit = Polynomial::Builder(table).degree(1).build();

        polyfit.save(fileName);
        Polynomial loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    SECTION("Polynomial of degree 4") {
        Polynomial polyfit = Polynomial::Builder(table).degree(4).build();

        polyfit.save(fileName);
        Polynomial loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    if (dim > 1) {
        SECTION("Polynomial of differing degrees in each dimension")
        {
            auto degrees = std::vector<unsigned int>(dim);
            for (unsigned int i = 0; i < dim; i++)
                degrees.at(i) = i + 1;

            Polynomial polyfit = Polynomial::Builder(table).degree(degrees).build();

            polyfit.save(fileName);
            Polynomial loadedPolyfit(fileName);

            REQUIRE(polyfit == loadedPolyfit);
        }
    }

    remove(fileName);
}

TEST_CASE("Polynomial can be saved and loaded (2)", COMMON_TAGS)
{
    unsigned int dim = 3;
    std::vector<unsigned int> degrees = {1,2,3};
    DenseVector coeffs = DenseVector::Ones(24);
    coeffs(5) = 6; coeffs(10) = 11; coeffs(15) = 16; coeffs(20) = 21;

    const char *fileName = "test.polynomial";

    SECTION("Polynomial of degrees (1,2,3)")
    {
        Polynomial poly(degrees, coeffs);
        poly.save(fileName);
        Polynomial loadedPoly(fileName);

        REQUIRE(poly == loadedPoly);
    }

    remove(fileName);
}

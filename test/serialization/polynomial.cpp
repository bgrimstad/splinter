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
        DenseMatrix powers(3, dim);
        powers << 1, 0,
            0, 1,
            1, 1;
        Polynomial polyfit = Polynomial::Builder(table).powers(powers).build();

        polyfit.save(fileName);
        Polynomial loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    SECTION("Polynomial of degree 2") {
        DenseMatrix powers(4, dim);
        powers << 1, 0,
            0, 1,
            0, 3,
            2, 0;
        Polynomial polyfit = Polynomial::Builder(table).powers(powers).build();

        polyfit.save(fileName);
        Polynomial loadedPolyfit(fileName);

        REQUIRE(polyfit == loadedPolyfit);
    }

    if (dim > 1) {
        SECTION("Polynomial with many terms")
        {
            int maxPower = 6;
            DenseMatrix powers(maxPower*maxPower, dim);

            for (int i = 0; i < maxPower; ++i)
            {
                for (int j = 0; j < maxPower; ++j)
                {
                    powers(i*maxPower + j, 0) = i;
                    powers(i*maxPower + j, 1) = j;
                }
            }

            Polynomial polyfit = Polynomial::Builder(table).powers(powers).build();

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
    DenseMatrix powers(3, dim);
    powers << 1, 2, 3,
        0, 1, 3,
        3, 0, 0;
    DenseVector coeffs = DenseVector::Ones(24);
    coeffs(5) = 6; coeffs(10) = 11; coeffs(15) = 16; coeffs(20) = 21;

    const char *fileName = "test.polynomial";

    SECTION("Polynomial of powers (1,2,3, 0,1,3, 3,0,0)")
    {
        Polynomial poly(powers, coeffs);
        poly.save(fileName);
        Polynomial loadedPoly(fileName);

        REQUIRE(poly == loadedPoly);
    }

    remove(fileName);
}

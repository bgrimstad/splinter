/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <utils/test_utils.h>
#include <utils/test_function_utils.h>
#include <utilities.h>
#include <bspline_builders.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][bspline]"
#define COMMON_TEXT " value approximation test with polynomials"


// TESTING LINEAR B-SPLINES
TEST_CASE("Linear BSpline function" COMMON_TEXT " densely sampled", COMMON_TAGS "[bsplinetype::linear][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 1);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 1.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             5000,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Linear BSpline function" COMMON_TEXT " sampled with medium density", COMMON_TAGS "[bsplinetype::linear][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 1);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.25;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 1.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             500,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Linear BSpline function" COMMON_TEXT " sparsely sampled", COMMON_TAGS "[bsplinetype::linear][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 1);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.6;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 1.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             50,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Linear BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::linear][jacobian]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 1);
        return (Function*) new BSpline(bs);
    };

    double one_eps = 5e-5;
    double two_eps = 5e-5;
    double inf_eps = 5e-5;

    for (auto test_func : getPolynomialFunctions())
    {
        compareJacobianValue(test_func, bspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TESTING QUADRATIC B-SPLINES
TEST_CASE("Quadratic BSpline function" COMMON_TEXT " densely sampled", COMMON_TAGS "[bsplinetype::quadratic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 2);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 2.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             5000,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quadratic BSpline function" COMMON_TEXT " sampled with normal density", COMMON_TAGS "[bsplinetype::quadratic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 2);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 2.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             500,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quadratic BSpline function" COMMON_TEXT " sparsely sampled", COMMON_TAGS "[bsplinetype::quadratic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 2);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.7;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 2.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             50,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quadratic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quadratic][jacobian]")
{
    auto bspline_builder =  [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 2);
        return (Function*) new BSpline(bs);
    };

    double one_eps = 6e-5;
    double two_eps = 6e-5;
    double inf_eps = 6e-5;

    for (auto test_func : getPolynomialFunctions())
    {
        compareJacobianValue(test_func, bspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TESTING CUBIC B-SPLINES
TEST_CASE("Cubic BSpline function" COMMON_TEXT " densely sampled", COMMON_TAGS "[bsplinetype::cubic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 3);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 3.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             5000,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Cubic BSpline function" COMMON_TEXT " sampled with normal density", COMMON_TAGS "[bsplinetype::cubic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 3);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 3.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             500,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Cubic BSpline function" COMMON_TEXT " sparsely sampled", COMMON_TAGS "[bsplinetype::cubic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 3);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.2;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 3.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             80,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Cubic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::cubic][jacobian]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 3);
        return (Function*) new BSpline(bs);
    };

    double one_eps = 6e-5;
    double two_eps = 6e-5;
    double inf_eps = 6e-5;

    for (auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc, bspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TESTING QUARTIC B-SPLINES
TEST_CASE("Quartic BSpline function" COMMON_TEXT " densely sampled", COMMON_TAGS "[bsplinetype::quartic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 4);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 4.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             5000,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quartic BSpline function" COMMON_TEXT " sampled with normal density", COMMON_TAGS "[bsplinetype::quartic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 4);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 4.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             500,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quartic BSpline function" COMMON_TEXT " sparsely sampled", COMMON_TAGS "[bsplinetype::quartic][function-value]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 4);
        return (Function*) new BSpline(bs);
    };

    for (auto test_func : getPolynomialFunctions())
    {
        double one_eps = 0.1;
        double two_eps = 0.1;
        double inf_eps = 0.1;

        // If the degree of the exact function is less than or equal to the degree
        // of the B-Spline we are using to approximate it, the B-Spline should approximate
        // the function exactly.
        if (test_func->isConstDegree() && test_func->getMaxDegree() <= 4.0)
        {
            one_eps = 1e-5;
            two_eps = 1e-5;
            inf_eps = 1e-5;
        }

        compareFunctionValue(test_func, bspline_builder,
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Quartic BSpline jacobian" COMMON_TEXT, COMMON_TAGS "[bsplinetype::quartic][jacobian]")
{
    auto bspline_builder = [](const DataTable &table)
    {
        BSpline bs = bspline_interpolator(table, 4);
        return (Function*) new BSpline(bs);
    };

    double one_eps = 7e-5;
    double two_eps = 7e-5;
    double inf_eps = 7e-5;

    for (auto test_func : getPolynomialFunctions())
    {
        compareJacobianValue(test_func, bspline_builder,
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}


// TESTING INTERPOLATION OF ZERO FUNCTION (SPARSE SOLVER)
TEST_CASE("B-spline approximation of zero function (sparse)", COMMON_TAGS "[zero-function][sparse]")
{
    DataTable data;

    auto x0 = linspace(-10, 10, 20);
    auto x1 = linspace(-10, 10, 20);

    for (auto x0_i : x0)
        for (auto x1_i : x1)
            data.add_sample({x0_i, x1_i}, .0);

    REQUIRE_NOTHROW(bspline_interpolator(data, 0));
    REQUIRE_NOTHROW(bspline_interpolator(data, 1));
    REQUIRE_NOTHROW(bspline_interpolator(data, 2));
    REQUIRE_NOTHROW(bspline_interpolator(data, 3));
    REQUIRE_NOTHROW(bspline_interpolator(data, 4));
    REQUIRE_NOTHROW(bspline_interpolator(data, 5));
}

// TESTING INTERPOLATION OF ZERO FUNCTION (DENSE SOLVER)
TEST_CASE("B-spline approximation of zero function (dense)", COMMON_TAGS "[zero-function][dense]")
{
    DataTable data;

    auto x0 = linspace(-10, 10, 8);
    auto x1 = linspace(-10, 10, 8);

    for (auto x0_i : x0)
        for (auto x1_i : x1)
            data.add_sample({x0_i, x1_i}, .0);

    REQUIRE_NOTHROW(bspline_interpolator(data, 0));
    REQUIRE_NOTHROW(bspline_interpolator(data, 1));
    REQUIRE_NOTHROW(bspline_interpolator(data, 2));
    REQUIRE_NOTHROW(bspline_interpolator(data, 3));
    REQUIRE_NOTHROW(bspline_interpolator(data, 4));
    REQUIRE_NOTHROW(bspline_interpolator(data, 5));
}
/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "term.h"
#include <Catch.h>
#include <testingutilities.h>
#include <radialbasisfunction.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomial][rbf][radialbasisfunction]"
#define COMMON_TEXT " value approximation test with polynomials"


/*
 * Gaussian
 */
TEST_CASE("Gaussian RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
    },
    TestType::FunctionValue,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Gaussian RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Gaussian RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}


/*
 * Inverse multiquadric
 */
TEST_CASE("Inverse Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
    },
    TestType::FunctionValue,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Inverse Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Inverse Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}


/*
 * Inverse quadric
 */
TEST_CASE("Inverse Quadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
    },
    TestType::FunctionValue,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Inverse Quadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Inverse Quadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}


/*
 * Multiquadric
 */
TEST_CASE("Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
    },
    TestType::FunctionValue,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}


/*
 * Thin plate spline
 */
TEST_CASE("Thin plate spline RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
    },
    TestType::FunctionValue,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Thin plate spline RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
    },
    TestType::Jacobian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

TEST_CASE("Thin plate spline RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    testApproximation(getPolynomialFunctions(),
        [](const DataTable &table) {
        return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
    },
    TestType::Hessian,
    300,  // Number of points to sample at
    1337, // Number of points to test against
    one_eps, two_eps, inf_eps);
}

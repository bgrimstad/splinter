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
#include <polynomialregression.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomialregression][polynomial]"
#define COMMON_TEXT " value approximation test with polynomials"


TEST_CASE("PolynomialRegression function" COMMON_TEXT, COMMON_TAGS "[function-value]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        INFO("PolynomialRegression of degree " << degree);
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::FunctionValue,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }

}

TEST_CASE("PolynomialRegression jacobian" COMMON_TEXT, COMMON_TAGS "[jacobian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        INFO("PolynomialRegression of degree " << degree);
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::Jacobian,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }

}

TEST_CASE("PolynomialRegression hessian" COMMON_TEXT, COMMON_TAGS "[hessian]") {
    // TODO: These should probably be global?
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(int degree = 1; degree < 10; ++degree) {
        INFO("PolynomialRegression of degree " << degree);
        testApproximation(getPolynomialFunctions(),
        [degree](const DataTable &table) {
            return (Approximant *) new PolynomialRegression(table, degree);
        },
        TestType::Hessian,
        300,  // Number of points to sample at
        1337, // Number of points to test against
        one_eps, two_eps, inf_eps);
    }

}

//TEST_CASE("PolynomialRegression of degree 1 function value approximation test with polynomials", "[approximation][polynomialregression][bsplinetype::linear][function-value][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::LINEAR);
//    },
//    TestType::FunctionValue,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Linear BSpline jacobian value approximation test with polynomials", "[approximation][bspline][bsplinetype::linear][jacobian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::LINEAR);
//    },
//    TestType::Jacobian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Linear BSpline hessian value approximation test with polynomials", "[approximation][bspline][bsplinetype::linear][hessian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::LINEAR);
//    },
//    TestType::Hessian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}



//TEST_CASE("Quadratic BSpline function value approximation test with polynomials", "[approximation][bspline][bsplinetype::quadratic][function-value][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
//    },
//    TestType::FunctionValue,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Quadratic BSpline jacobian value approximation test with polynomials", "[approximation][bspline][bsplinetype::quadratic][jacobian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
//    },
//    TestType::Jacobian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Quadratic BSpline hessian value approximation test with polynomials", "[approximation][bspline][bsplinetype::quadratic][hessian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUADRATIC);
//    },
//    TestType::Hessian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}



//TEST_CASE("Cubic BSpline function value approximation test with polynomials", "[approximation][bspline][bsplinetype::cubic][function-value][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::CUBIC);
//    },
//    TestType::FunctionValue,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Cubic BSpline jacobian value approximation test with polynomials", "[approximation][bspline][bsplinetype::cubic][jacobian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::CUBIC);
//    },
//    TestType::Jacobian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Cubic BSpline hessian value approximation test with polynomials", "[approximation][bspline][bsplinetype::cubic][hessian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::CUBIC);
//    },
//    TestType::Hessian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}



//TEST_CASE("Quartic BSpline function value approximation test with polynomials", "[approximation][bspline][bsplinetype::quartic][function-value][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
//    },
//    TestType::FunctionValue,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Quartic BSpline jacobian value approximation test with polynomials", "[approximation][bspline][bsplinetype::quartic][jacobian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
//    },
//    TestType::Jacobian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}

//TEST_CASE("Quartic BSpline hessian value approximation test with polynomials", "[approximation][bspline][bsplinetype::quartic][hessian][polynomial]") {
//    // TODO: These should probably be global?
//    double one_eps = 0.1;
//    double two_eps = 0.1;
//    double inf_eps = 0.1;

//    testApproximation(getPolynomialFunctions(),
//        [](const DataTable &table) {
//        return (Approximant *) new BSpline(table, BSplineType::QUARTIC);
//    },
//    TestType::Hessian,
//    300,  // Number of points to sample at
//    1337, // Number of points to test against
//    one_eps, two_eps, inf_eps);
//}


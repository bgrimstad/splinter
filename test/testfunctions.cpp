/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "testfunctions.h"
#include "testfunction.h"

#define MAKE_TESTFUNCTION(FUNC, NUMVARIABLES, PRETTY_STRING) new TestFunction([](const std::vector<double> &x) { return FUNC; }, NUMVARIABLES, PRETTY_STRING)
#define MAKE_CONSTDEGREE_TESTFUNCTION(FUNC, NUMVARIABLES, PRETTY_STRING, DEGREES) new TestFunction([](const std::vector<double> &x) { return FUNC; }, NUMVARIABLES, PRETTY_STRING, DEGREES)

using namespace SPLINTER;


std::vector<std::vector<TestFunction *>> testFunctions = std::vector<std::vector<TestFunction *>>();


void setupTestFunctions() {
    // f_x_y: function of x variables and y degrees

    // Functions of one variable
    std::vector<double> one_degrees(1);
    one_degrees = {0};
    auto f_1_0 = MAKE_CONSTDEGREE_TESTFUNCTION(-13.37, 1, "-13.37", one_degrees);
    one_degrees = {1};
    auto f_1_1 = MAKE_CONSTDEGREE_TESTFUNCTION(-5.1*x.at(0) + 13.37, 1, "-5.1*x + 13.37", one_degrees);
    one_degrees = {2};
    auto f_1_2 = MAKE_CONSTDEGREE_TESTFUNCTION(8.1*x.at(0)*x.at(0) - 0.2*x.at(0) + 2313.1, 1, "8.1*x^2 - 0.2*x + 2313.1", one_degrees);
    one_degrees = {3};
    auto f_1_3 = MAKE_CONSTDEGREE_TESTFUNCTION(-4.5*x.at(0)*x.at(0)*x.at(0) + 2.2*x.at(0)*x.at(0), 1, "-4.5*x^3 + 2.2*x^2", one_degrees);
    one_degrees = {4};
    auto f_1_4 = MAKE_CONSTDEGREE_TESTFUNCTION(4.5*x.at(0)*x.at(0)*x.at(0)*x.at(0) + 3*x.at(0)*x.at(0)*x.at(0) - x.at(0), 1, "4.5*x^4 + x^3 - x^2", one_degrees);

    // Functions of two variables
    std::vector<double> two_degrees(2);
    two_degrees = {0,0};
    auto f_2_0 = MAKE_CONSTDEGREE_TESTFUNCTION(0.1, 2, "0.1", two_degrees);
    two_degrees = {1,1};
    auto f_2_1 = MAKE_CONSTDEGREE_TESTFUNCTION(- 5.1*x.at(0) + 13.37*x.at(1), 2, "-5.1*x + 13.37*y", two_degrees);
    two_degrees = {2,2};
    auto f_2_2 = MAKE_CONSTDEGREE_TESTFUNCTION(8.1*x.at(0)*x.at(0) - 0.2*x.at(0)*x.at(1) + 13.37*x.at(1)*x.at(1), 2, "8.1*x^2 - 0.2*x*y + 13.37*y^2", two_degrees);
    two_degrees = {3,2};
    auto f_2_3 = MAKE_CONSTDEGREE_TESTFUNCTION(-4.5*x.at(0)*x.at(0)*x.at(0) + 2.2*x.at(0)*x.at(0) - x.at(1)*x.at(1) + 3, 2, "-4.5*x^3 + 2.2*x^2 - y^2 + 3", two_degrees);
    two_degrees = {4,1};
    auto f_2_4 = MAKE_CONSTDEGREE_TESTFUNCTION(4.5*x.at(0)*x.at(0)*x.at(0)*x.at(0) - x.at(0)*x.at(0)*x.at(0) + 3*x.at(0)*x.at(0)*x.at(1), 2, "4.5*x^4 - x^3 + 3*x^2*y", two_degrees);
    two_degrees = {3,2};
    auto f_2_5 = MAKE_CONSTDEGREE_TESTFUNCTION(-57*x.at(1)*x.at(1)*x.at(0)*x.at(0)*x.at(0) - 0.1*x.at(0)*x.at(0)*x.at(0)*x.at(1)*x.at(1) + 1.1*x.at(1)*x.at(0)*x.at(0) + x.at(1) - 1e6, 2, "-57*x^3*y^2 - 0.1*x^3*y^2 + 1.1*x^2*y + y - 1e6", two_degrees);
    two_degrees = {6,4};
    // Six-hump camel back function
    auto f_2_6 = MAKE_CONSTDEGREE_TESTFUNCTION((4 - 2.1*x.at(0)*x.at(0) + (1/3.)*x.at(0)*x.at(0)*x.at(0)*x.at(0))*x.at(0)*x.at(0) + x.at(0)*x.at(1) + (-4 + 4*x.at(1)*x.at(1))*x.at(1)*x.at(1), 2, "(4 - 2.1*x^2 + 1/3*x^4) * x^2 + x*y + (-4 + 4*y^2)*y^2", two_degrees);

    // Functions of three variables
    // Note: Remove the . in the constant and several tests fail
    std::vector<double> three_degrees(3);
    three_degrees = {0,0,0};
    auto f_3_0 = MAKE_CONSTDEGREE_TESTFUNCTION(6.534460297, 3, "6.534460297", three_degrees);
    three_degrees = {1,1,1};
    auto f_3_1 = MAKE_CONSTDEGREE_TESTFUNCTION(x.at(0)+x.at(1)-x.at(2)-1, 3, "x + y - z - 1", three_degrees);

    /* The next two functions causes a lot of trouble for some of the Approximants, especially RBF and PolynomialRegression
     Main problem is probably that it is a three-dimensional function, so the number of grid points per dimension
     becomes small.
     */
    three_degrees = {1,2,2};
    auto f_3_2 = MAKE_CONSTDEGREE_TESTFUNCTION(x.at(1)*x.at(2) + x.at(1)*x.at(1) + 17.0*x.at(2)*x.at(2) - x.at(0) - 10, 3, "y*z + y^2 + 17.0*z^2 - x - 10", three_degrees);
    three_degrees = {2,3,3};
    auto f_3_3 = MAKE_CONSTDEGREE_TESTFUNCTION((x.at(0)+x.at(1)-x.at(2)-1) * (x.at(1)*x.at(2) + x.at(1)*x.at(1) + 17.0*x.at(2)*x.at(2) - x.at(0) - 10), 3, "(x + y - z - 1) * (y*z + y^2 + 17.0*z^2 - x - 10)", three_degrees);

    // Non-polynomial (aka. nasty) functions
    auto f_nasty_0 = MAKE_TESTFUNCTION(std::exp(x.at(2)) * (std::pow(1.3, x.at(0)) + std::pow(x.at(1), x.at(0))), 3, "e^z * (1.3^x + y^x)");
    auto f_nasty_1 = MAKE_TESTFUNCTION(3*x.at(0)/(x.at(0)*x.at(0) + 3*x.at(1)*x.at(1)*x.at(1)), 2, "3*x / (x^2 + 3*y^3)");
    auto f_nasty_2 = MAKE_TESTFUNCTION(std::sin(x.at(0))*std::exp(x.at(0))*std::log(x.at(1))+std::cos(x.at(1))*std::tan(x.at(0)), 2, "sin(x)*e^x*log(y) + cos(y)*tan(x)");

    // First vector is the vector of nasty functions
    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(0).push_back(f_nasty_0);
    testFunctions.at(0).push_back(f_nasty_1);
    testFunctions.at(0).push_back(f_nasty_2);

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(1).push_back(f_1_0);
    testFunctions.at(1).push_back(f_1_1);
    testFunctions.at(1).push_back(f_1_2);
    testFunctions.at(1).push_back(f_1_3);
    testFunctions.at(1).push_back(f_1_4);

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(2).push_back(f_2_0);
    testFunctions.at(2).push_back(f_2_1);
    testFunctions.at(2).push_back(f_2_2);
    testFunctions.at(2).push_back(f_2_3);
    testFunctions.at(2).push_back(f_2_4);
    testFunctions.at(2).push_back(f_2_5);
    testFunctions.at(2).push_back(f_2_6);

    testFunctions.push_back(std::vector<TestFunction *>());
    testFunctions.at(3).push_back(f_3_0);
    testFunctions.at(3).push_back(f_3_1);
    testFunctions.at(3).push_back(f_3_2);
    testFunctions.at(3).push_back(f_3_3);
}

void tearDownTestFunctions() {
    for(auto &vec : testFunctions) {
        for(auto &testFunc : vec) {
            delete testFunc;
        }
    }
}

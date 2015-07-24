/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "testfunctions.h"
#include <term.h>
#include <termfunction.h>


std::vector<std::vector<TermFunction *>> testFunctions = std::vector<std::vector<TermFunction *>>();


void setupTestFunctions() {
    // f_x_y: function of x variables and y degrees

    Var x(0, "x");
    Var y(1, "y");
    Var z(2, "z");
    Var x1(3, "x1");
    Var y1(4, "y1");
    Var z1(5, "z1");
    Var x2(6, "x2");
    Var y2(7, "y2");
    Var z2(8, "z2");

    // Functions of one variable
    auto f_1_0 = -13.37 + 0*x;
    auto f_1_1 = -5.1*x + 13.37;
    auto f_1_2 = 8.1*(x^2) - 0.2*x + 2313.1;
    auto f_1_3 = -4.5*(x^3) + 2.2*(x^2);
    auto f_1_4 = x*(4.5*(x^3) + 3*(x^2) - x);

    // Functions of two variables
    auto f_2_0 = 0.1 + 0*x*y;
    auto f_2_1 = - 5.1*x + 13.37*y;
    auto f_2_2 = 8.1*(x^2) - 0.2*x*y + 13.37*(y^2);
    auto f_2_3 = -4.5*(x^3) + 2.2*(x^2) - (y^2) + 3;
    auto f_2_4 = x*(4.5*(x^3) - (x^2) + 3*x*y);
    auto f_2_5 = -57*(y^2)*(x^3) - 0.1*(x^3)*(y^2) + 1.1*y*(x^2) + y - 1e10;
    // Six-hump camel back function
    auto f_2_6 = (4 - 2.1*(x^2) + (1/3.)*(x^4))*(x^2) + x*y + (-4 + 4*(y^2))*(y^2);

    // Functions of three variables
    auto f_3_0 = 6534460297 + 0*x*y*z;
    auto f_3_1 = x+y-z-1;
    auto f_3_2 = y*z + 3.9*(z^2) + (y^2) + 13.1*(z^2) - x - 10;
    auto f_3_3 = f_3_2 * f_3_1;

    // Non-polynomial (aka. nasty) functions
    auto f_nasty_0 = E(z) * ((1.3^x) + (y^x));
    auto f_nasty_1 = 3*x/((x^2) + 3*(y^3));
    auto f_nasty_2 = Sin(x)*E(x)*Log(y)+Cos(y)*Tan(x);

    // First vector is the vector of nasty functions
    testFunctions.push_back(std::vector<TermFunction *>());
    testFunctions.at(0).push_back(new TermFunction(f_nasty_0));
    testFunctions.at(0).push_back(new TermFunction(f_nasty_1));
    testFunctions.at(0).push_back(new TermFunction(f_nasty_2));

    testFunctions.push_back(std::vector<TermFunction *>());
    testFunctions.at(1).push_back(new TermFunction(f_1_0));
    testFunctions.at(1).push_back(new TermFunction(f_1_1));
    testFunctions.at(1).push_back(new TermFunction(f_1_2));
    testFunctions.at(1).push_back(new TermFunction(f_1_3));
    testFunctions.at(1).push_back(new TermFunction(f_1_4));

    testFunctions.push_back(std::vector<TermFunction *>());
    testFunctions.at(2).push_back(new TermFunction(f_2_0));
    testFunctions.at(2).push_back(new TermFunction(f_2_1));
    testFunctions.at(2).push_back(new TermFunction(f_2_2));
    testFunctions.at(2).push_back(new TermFunction(f_2_3));
    testFunctions.at(2).push_back(new TermFunction(f_2_4));
    testFunctions.at(2).push_back(new TermFunction(f_2_5));
    testFunctions.at(2).push_back(new TermFunction(f_2_6));

    testFunctions.push_back(std::vector<TermFunction *>());
    testFunctions.at(3).push_back(new TermFunction(f_3_0));
    testFunctions.at(3).push_back(new TermFunction(f_3_1));
    testFunctions.at(3).push_back(new TermFunction(f_3_2));
    testFunctions.at(3).push_back(new TermFunction(f_3_3));
}

void tearDownTestFunctions() {
    for(auto &vec : testFunctions) {
        for(auto &term : vec) {
            delete term;
        }
    }
}

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <testfunction.h>
#include <test_functions.h>
#include <testingutilities.h>

using namespace SPLINTER;

TestFunction::TestFunction(unsigned int n, Term *func)
        : numVariables(n),
          f(func->simplify())
{
    std::cout << "f: ";
    f->pretty_text(std::cout); std::cout << std::endl;
    calculateJacobian();
    calculateHessian();
}
TestFunction::TestFunction(unsigned int n, Term &func)
        : TestFunction(n, func.clone())
{
}

double TestFunction::eval(DenseVector x) const
{
    auto xvec = denseToVec(x);

    return f->eval(xvec);
}

DenseMatrix TestFunction::evalJacobian(DenseVector x) const
{
    DenseMatrix jacobian(1, numVariables);

    auto xvec = denseToVec(x);

    for(size_t i = 0; i < numVariables; i++) {
        jacobian(0, i) = jac.at(i)->eval(xvec);
    }

    return jacobian;
}

DenseMatrix TestFunction::evalHessian(DenseVector x) const
{
    DenseMatrix hessian(numVariables, numVariables);

    auto xvec = denseToVec(x);

    for(size_t i = 0; i < numVariables; i++) {
        for(size_t j = 0; j < numVariables; j++) {
            hessian(i, j) = hes.at(i).at(j)->eval(xvec);
        }
    }

    return hessian;
}

void TestFunction::calculateJacobian()
{
    jac.clear();
    for(size_t i = 0; i < numVariables; i++) {
        auto df = f->derivative(i);

        std::cout << "df/dx: ";
        df->pretty_text(std::cout);
        std::cout << std::endl;

        auto simplified_df = df->simplify();
        std::cout << "simplified df/dx" << i << ": ";
        simplified_df->pretty_text(std::cout);
        std::cout << std::endl;
        std::cout << std::endl;

        jac.push_back(simplified_df);
    }
}

void TestFunction::calculateHessian()
{
    hes.clear();

    for(size_t i = 0; i < numVariables; i++) {
        hes.push_back(std::vector<Term *>());

        for(size_t j = 0; j < numVariables; j++) {
            auto ddf = jac.at(i)->derivative(j);

            std::cout << "d^2f/dx" << i << "dx" << j << ": ";
            ddf->pretty_text(std::cout);
            std::cout << std::endl;

            auto simplified_ddf = ddf->simplify();
            std::cout << "simplified d^2f/dx" << i << "dx" << j << ": ";
            simplified_ddf->pretty_text(std::cout);
            std::cout << std::endl;
            std::cout << std::endl;

            hes.at(i).push_back(simplified_ddf);
        }
    }
}

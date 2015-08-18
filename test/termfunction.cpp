/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "termfunction.h"
#include "term.h"
#include <testingutilities.h>

using namespace SPLINTER;

TermFunction::TermFunction(Term *func)
        : f(func)
{
    // Gather variables before simplification so we don't lose the number
    // of variables. If a func: f = 0*x + y + 13 is passed in,
    // we still want it to be treated as a function of two variables, not one.
    gatherVariables();

    auto temp = f->simplify();
    if(temp != f) {
        delete f;
        f = temp;
    }

    calculateJacobian();
    calculateHessian();

#ifndef NDEBUG
    printAll();
#endif // ifndef NDEBUG
}

TermFunction::TermFunction(Term &func)
        : TermFunction(func.clone())
{
}

TermFunction::TermFunction(Term &&func)
        : TermFunction(func.clone())
{
}

TermFunction::~TermFunction()
{
    delete f;
    for(auto &df : jac) {
        delete df;
    }
    for(auto &ddfvec : hes) {
        for(auto &ddf : ddfvec) {
            delete ddf;
        }
    }
}

void TermFunction::gatherVariables()
{
    std::set<Var> temp;
    f->getVariables(temp);

    for(auto &var : temp) {
        variables.push_back(var);
    }

    sort(variables.begin(), variables.end());
    numVariables = variables.size();

#ifndef NDEBUG
    std::cout << "NumVariables: " << numVariables << std::endl;
#endif
}

double TermFunction::eval(DenseVector x) const
{
    auto xvec = denseToVec(x);

    return f->eval(xvec);
}

DenseMatrix TermFunction::evalJacobian(DenseVector x) const
{
    DenseMatrix jacobian(1, numVariables);

    auto xvec = denseToVec(x);

    for(size_t i = 0; i < numVariables; i++) {
        jacobian(i) = jac.at(i)->eval(xvec);
    }

    return jacobian;
}

DenseMatrix TermFunction::evalHessian(DenseVector x) const
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

void TermFunction::calculateJacobian()
{
    jac.clear();
    for(size_t i = 0; i < numVariables; i++) {
        auto dx = variables.at(i);
        jac.push_back(f->derivative(dx));
    }
}

// TODO: Take advantage of the fact that the Hessian is symmetric
// over the diagonal. If we do this we can skip almost half the differentations
void TermFunction::calculateHessian()
{
    hes.clear();

    for(size_t i = 0; i < numVariables; i++) {
        hes.push_back(std::vector<Term *>());

        for(size_t j = 0; j < numVariables; j++) {
            Var dx = variables.at(j);
            hes.at(i).push_back(jac.at(i)->derivative(dx));
        }
    }
}

void TermFunction::printAll() const
{
    std::cout << getFunctionStr();
    std::cout << getJacobianStr();
    std::cout << getHessianStr();
    std::cout << std::endl;
}

std::string TermFunction::getFunctionStr() const
{
    std::stringstream s;
    s << "Function: ";
    f->pretty_text(s);
    s << std::endl;
    return s.str();
}

std::string TermFunction::getJacobianStr() const
{
    std::stringstream s;
    s << "Jacobian: " << std::endl;
    for(size_t i = 0; i < numVariables; i++) {
        s << "df/d";
        variables.at(i).pretty_text(s);
        s << ": ";
        jac.at(i)->pretty_text(s);
        s << std::endl;
    }
    return s.str();
}

std::string TermFunction::getHessianStr() const
{
    std::stringstream s;
    s << "Hessian: " << std::endl;
    for(size_t i = 0; i < numVariables; i++) {
        for(size_t j = 0; j < numVariables; j++) {
            s << "d^2f/d";
            variables.at(i).pretty_text(s);
            if(i == j) {
                s << "^2";
            } else {
                s << "d";
                variables.at(j).pretty_text(s);
            }
            s << ": ";
            hes.at(i).at(j)->pretty_text(s);
            s << std::endl;
        }
    }
    return s.str();
}
/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_TESTFUNCTION_H
#define SPLINTER_TESTFUNCTION_H

#include <function.h>

using namespace SPLINTER;

class Term;
class Var;

// TODO: Keep track of the domain of the function
// This is very useful for testing functions that have "interesting" regions
class TermFunction : public Function
{
public:
    TermFunction(Term *func);
    TermFunction(Term &func);
    TermFunction(Term &&func);

    virtual ~TermFunction();

    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override;

    inline Term *getF() const { return f; }
    inline std::vector<Term *> getDF() const { return jac; }
    inline std::vector<std::vector<Term *>> getDDF() const { return hes; }

    inline unsigned int getNumVariables() const override { return numVariables; }

    void printAll() const;

private:
    unsigned int numVariables;
    std::vector<Var> variables;

    Term *f;
    std::vector<Term *> jac;
    std::vector<std::vector<Term *>> hes;

    void gatherVariables();
    void calculateJacobian();
    void calculateHessian();
};

#endif // SPLINTER_TESTFUNCTION_H

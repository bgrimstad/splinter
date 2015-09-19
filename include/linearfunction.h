/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_LINEARFUNCTION_H
#define SPLINTER_LINEARFUNCTION_H

#include "definitions.h"
#include "function.h"

namespace SPLINTER
{

/**
 * Interface for functions (linear in the coefficients) on the form:
 * f(x) = c(1)*b1(x) + c(2)*b2(x) + ... + c(n)*bn(x) = c^T * b(x),
 * where c are coefficients and b is a vector of basis functions.
 */
class SPLINTER_API LinearFunction : public Function
{
public:
    virtual ~LinearFunction() {}

    //virtual DenseVector computeCoefficients(DataTable data);
    unsigned int getNumCoefficients() const { return coefficients.size(); }
protected:
    LinearFunction(unsigned int numVariables)
            : Function(numVariables) {}

    DenseVector coefficients;
};

} // namespace SPLINTER

#endif // SPLINTER_LINEARFUNCTION_H

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
 * Interface for functions (linear in the coefficients):
 * f(x) = c(1)*b1(x) + c(2)*b2(x) + ... + c(n)*bn(x) = c^T * b(x),
 * where c are coefficients and b is a vector of basis functions.
 */
template <class Vec, class Mat>
class SPLINTER_API LinearFunction : public Function
{
public:
    virtual ~LinearFunction() {}

    /**
     * Returns the function value at x
     */
    double eval(DenseVector x) const override
    {
        checkInput(x);
        auto res = coefficients.transpose()*evalBasis(x);
        return res(0);
    }

    /**
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    DenseMatrix evalJacobian(DenseVector x) const override
    {
        checkInput(x);
        return coefficients.transpose()*evalBasisJacobian(x);
    }

    /**
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     * TODO: implement
     */
    DenseMatrix evalHessian(DenseVector x) const override
    {
        throw Exception("LinearFunction::evalHessian: Not implemented!");
//        return DenseMatrix::Zero(numVariables, numVariables);
    }

    /**
     * Evaluate basis functions at x
     */
    virtual Vec evalBasis(DenseVector x) const = 0;

    /**
     * Evaluate Jacobian of basis functions at x
     */
    virtual Mat evalBasisJacobian(DenseVector x) const = 0;

    DenseVector getCoefficients()
    {
        return coefficients;
    }

    void setCoefficients(DenseVector coefficients)
    {
        this->coefficients = coefficients;
    }

    unsigned int getNumCoefficients() const {
        return coefficients.size();
    }

protected:
    LinearFunction(unsigned int numVariables, DenseVector coefficients)
            : Function(numVariables),
              coefficients(coefficients)
    {}

    DenseVector coefficients;

    friend class Serializer;

//    friend bool operator== <>(const LinearFunction<Vec, Mat> &lhs, const LinearFunction<Vec, Mat> &rhs);
};


} // namespace SPLINTER

#endif // SPLINTER_LINEARFUNCTION_H

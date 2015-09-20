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

    /**
     * Returns the function value at x
     */
    double eval(DenseVector x) const override;

    /**
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    DenseMatrix evalJacobian(DenseVector x) const override;

    /**
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     * TODO: implement
     */
    DenseMatrix evalHessian(DenseVector x) const override
    {
        return DenseMatrix::Zero(numVariables, numVariables);
    }

    /**
     * Evaluate basis functions at x
     */
    virtual SparseVector evalBasisFunctions(DenseVector x) const = 0;

    /**
     * Evaluate Jacobian of basis functions at x
     */
    virtual SparseMatrix evalBasisFunctionsJacobian(DenseVector x) const = 0;

    DenseMatrix getCoefficients()
    {
        return coefficients;
    }

    void setCoefficients(DenseMatrix coefficients)
    {
        if (coefficients.cols() != 1)
            throw Exception("LinearFunction::setCoefficients: coefficient matrix must have exactly one column.");
        this->coefficients = coefficients;
    }

    unsigned int getNumCoefficients() const { return coefficients.size(); }

protected:
    LinearFunction(unsigned int numVariables, DenseMatrix coefficients)
            : Function(numVariables),
              coefficients(coefficients)
    {}

    /**
     * Coefficients
     * NOTE: serialization fails when using DenseVector
     */
    DenseMatrix coefficients;

    friend class Serializer;
    friend bool operator==(const LinearFunction &lhs, const LinearFunction &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_LINEARFUNCTION_H

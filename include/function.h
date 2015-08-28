/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_FUNCTION_H
#define SPLINTER_FUNCTION_H

#include "definitions.h"

namespace SPLINTER
{

/*
 * Interface for functions
 */
class SPLINTER_API Function
{
public:
    Function()
        : Function(1) {}

    Function(unsigned int numVariables)
        : numVariables(numVariables) {}

    virtual ~Function() {}

    /**
     * Returns the function value at x
     */
    virtual double eval(DenseVector x) const;

    /**
     * Returns the function value at x
     */
    virtual double eval(const std::vector<double> &x) const;

    /**
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    virtual DenseMatrix evalJacobian(DenseVector x) const;

    /**
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    virtual std::vector<double> evalJacobian(const std::vector<double> &x) const;

    /**
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     */
    virtual DenseMatrix evalHessian(DenseVector x) const;

    /**
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     */
    virtual std::vector<std::vector<double>> evalHessian(const std::vector<double> &x) const;

    /**
     * Get the dimension
     */
    inline unsigned int getNumVariables() const
    {
        return numVariables;
    }

    /**
     * Returns the central difference at x
     * Vector of numVariables length
     */
    std::vector<double> centralDifference(const std::vector<double> &x) const;
    DenseMatrix centralDifference(DenseVector x) const;

    std::vector<std::vector<double>> secondOrderCentralDifference(const std::vector<double> &x) const;
    DenseMatrix secondOrderCentralDifference(DenseVector x) const;

protected:
    unsigned int numVariables; // Dimension of domain (size of x)

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_FUNCTION_H
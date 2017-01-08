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
#include "saveable.h"

namespace SPLINTER
{

/*
 * Interface for functions f : R^m -> R^n
 * All functions working with standard C++11 types are defined in terms of their Eigen counterparts.
 * Default implementations of Jacobian and Hessian evaluation is using central difference.
 * TODO: Remove current requirement that all functions must implement save and load!
 */
class SPLINTER_API Function : public Saveable
{
public:
    Function()
        : Function(1, 1) {}

    Function(unsigned int m, unsigned int n)
        : dimX(m), dimY(n) {}

    virtual ~Function() {}

    /**
     * Returns the function value at x
     */
    virtual std::vector<double> eval(const std::vector<double> &x) const = 0;

    virtual DenseVector eval(const DenseVector &x) const;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    virtual DenseMatrix evalJacobian(const DenseVector &x) const;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    std::vector<std::vector<double>> evalJacobian(const std::vector<double> &x) const;

    /**
     * Get dimensions
     */
    inline unsigned int getDimX() const
    {
        return dimX;
    }

    inline unsigned int getDimY() const
    {
        return dimY;
    }

    /**
     * Check input
     */
    void checkInput(const std::vector<double> &x) const {
        if (x.size() != dimX)
            throw Exception("Function::checkInput: Wrong dimension on evaluation point x.");
    }

    void checkInput(const DenseVector &x) const;

    /**
     * Returns the central difference at x
     */
    DenseMatrix centralDifference(const DenseVector &x) const;

    /**
     * Description of function.
     */
    virtual std::string getDescription() const
    {
        return "";
    }

protected:
    unsigned int dimX; // Dimension of domain (size of x)
    unsigned int dimY; // Dimension of codomain (size of y)

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_FUNCTION_H
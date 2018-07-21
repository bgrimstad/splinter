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
 * Interface for functions f : R^m -> R^n
 * All functions working with standard C++11 types are defined in terms of their Eigen counterparts.
 * Default implementations of Jacobian uses central differences.
 */
class SPLINTER_API Function
{
public:
    Function()
        : Function(1, 1) {}

    Function(unsigned int m, unsigned int n)
        : dim_x(m), dim_y(n) {}

    virtual ~Function() {}

    /**
     * Returns the (dimY) function values at x
     */
    virtual std::vector<double> eval(const std::vector<double> &x) const = 0;

    virtual DenseVector eval(const DenseVector &x) const;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    virtual DenseMatrix eval_jacobian(const DenseVector &x) const;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    std::vector<std::vector<double>> eval_jacobian(const std::vector<double> &x) const;

    /**
     * Get dimensions
     */
    inline unsigned int get_dim_x() const
    {
        return dim_x;
    }

    inline unsigned int get_dim_y() const
    {
        return dim_y;
    }

    /**
     * Check input
     */
    void check_input(const std::vector<double> &x) const {
        if (x.size() != dim_x)
            throw Exception("Function::check_input: Wrong dimension on evaluation point x.");
    }

    void check_input(const DenseVector &x) const;

    /**
     * Returns the central difference at x
     */
    DenseMatrix central_difference(const DenseVector &x) const;

    /**
     * Description of function.
     */
    virtual std::string get_description() const
    {
        return "";
    }

protected:
    unsigned int dim_x; // Dimension of domain (size of x)
    unsigned int dim_y; // Dimension of codomain (size of y)

};

} // namespace SPLINTER

#endif // SPLINTER_FUNCTION_H
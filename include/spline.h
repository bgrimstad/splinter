/*
This file is part of the Multivariate Splines library.
Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#ifndef SPLINE_H
#define SPLINE_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

/*
 * Interface for spline classes
 */
class Spline
{
public:
    Spline() {}
    virtual ~Spline() {}

    /*
     * Returns the spline value at x
     */
    virtual double eval(DenseVector &x) const = 0;

    /*
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    virtual DenseMatrix evalJacobian(DenseVector &x) const = 0;

    /*
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     */
    virtual DenseMatrix evalHessian(DenseVector &x) const = 0;

};

} // namespace MultivariateSplines

#endif // SPLINE_H

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_SPLINE_H
#define SPLINTER_SPLINE_H

#include "generaldefinitions.h"

namespace SPLINTER
{

/*
 * Interface for approximants
 */
class Approximant
{
public:
    Approximant() {}
    virtual ~Approximant() {}

    /*
     * Returns the spline value at x
     */
    virtual double eval(DenseVector x) const = 0;

    /*
     * Returns the (1 x numVariables) Jacobian evaluated at x
     */
    virtual DenseMatrix evalJacobian(DenseVector x) const = 0;

    /*
     * Returns the (numVariables x numVariables) Hessian evaluated at x
     */
    virtual DenseMatrix evalHessian(DenseVector x) const = 0;

    /*
    * Serialize and save spline to fileName
    * Throws if file could not be opened
    */
    virtual void save(const std::string fileName) const = 0;

    /*
    * Deserialize and load spline from fileName
    * Throws if file could not be opened or if the file format is wrong
    */
    virtual void load(const std::string fileName) = 0;
};

} // namespace SPLINTER

#endif // SPLINTER_SPLINE_H

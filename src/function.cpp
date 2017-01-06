/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <function.h>
#include "utilities.h"

namespace SPLINTER
{

std::vector<double> Function::eval(const std::vector<double> &x) const {
    return eval(stdToEigVec(x));
}

std::vector<std::vector<double>> Function::evalJacobian(const std::vector<double> &x) const
{
    auto denseX = stdToEigVec(x);

    return eigMatToStdVecVec(evalJacobian(denseX));
}

DenseMatrix Function::evalJacobian(const DenseVector &x) const
{
    return centralDifference(x);
}

DenseMatrix Function::centralDifference(const DenseVector &x) const
{
    DenseMatrix Jac(dimY, dimX);

    double h = 1e-6; // perturbation step size
    double hForward = 0.5*h;
    double hBackward = 0.5*h;

    for (unsigned int i = 0; i < dimX; ++i)
    {
        DenseVector xForward(x);
        xForward(i) = xForward(i) + hForward;

        DenseVector xBackward(x);
        xBackward(i) = xBackward(i) - hBackward;

        auto yForward = eval(xForward);
        auto yBackward = eval(xBackward);

        for (unsigned int j = 0; j < dimY; ++j)
            Jac(j, i) = (yForward.at(j) - yBackward.at(j)) / (hBackward + hForward);
    }

    return Jac;
}

void Function::checkInput(const DenseVector &x) const {
    return checkInput(eigToStdVec(x));
}


} // namespace SPLINTER
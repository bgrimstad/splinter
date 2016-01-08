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

double Function::eval(const std::vector<double> &x) const
{
    auto denseX = vectorToDenseVector(x);

    return eval(denseX);
}

std::vector<double> Function::evalJacobian(const std::vector<double> &x) const
{
    return centralDifference(x);
}

/**
 * Returns the second order central difference in x
 */
std::vector<std::vector<double>> Function::evalHessian(const std::vector<double> &x) const
{
    return secondOrderCentralDifference(x);
}

std::vector<double> Function::centralDifference(const std::vector<double> &x) const
{
    std::vector<double> dx(x.size());

    double h = 1e-6; // perturbation step size
    double hForward = 0.5*h;
    double hBackward = 0.5*h;

    for (size_t i = 0; i < getNumVariables(); ++i)
    {
        std::vector<double> xForward(x);
//        if (xForward(i) + hForward > variables.at(i)->getUpperBound())
//        {
//            hForward = 0;
//        }
//        else
//        {
        xForward.at(i) = xForward.at(i) + hForward;
//        }

        std::vector<double> xBackward(x);
//        if (xBackward(i) - hBackward < variables.at(i)->getLowerBound())
//        {
//            hBackward = 0;
//        }
//        else
//        {
        xBackward.at(i) = xBackward.at(i) - hBackward;
//        }

        double yForward = eval(xForward);
        double yBackward = eval(xBackward);

        dx.at(i) = (yForward - yBackward)/(hBackward + hForward);
    }

    return dx;
}

std::vector<std::vector<double>> Function::secondOrderCentralDifference(const std::vector<double> &x) const
{
    std::vector<std::vector<double>> ddx(getNumVariables());

    double h = 1e-6; // perturbation step size
    double hForward = 0.5*h;
    double hBackward = 0.5*h;

    for (size_t i = 0; i < getNumVariables(); ++i)
    {
        ddx.at(i) = std::vector<double>(getNumVariables());

        for (size_t j = 0; j < getNumVariables(); ++j)
        {
            std::vector<double> x0(x);
            std::vector<double> x1(x);
            std::vector<double> x2(x);
            std::vector<double> x3(x);

            x0.at(i) = x0.at(i) + hForward;
            x0.at(j) = x0.at(j) + hForward;

            x1.at(i) = x1.at(i) - hBackward;
            x1.at(j) = x1.at(j) + hForward;

            x2.at(i) = x2.at(i) + hForward;
            x2.at(j) = x2.at(j) - hBackward;

            x3.at(i) = x3.at(i) - hBackward;
            x3.at(j) = x3.at(j) - hBackward;

            ddx.at(i).at(j) = (eval(x0) - eval(x1) - eval(x2) + eval(x3)) / (h * h);
        }
    }

    return ddx;
}

DenseMatrix Function::secondOrderCentralDifference(DenseVector x) const
{
    auto vec = denseVectorToVector(x);

    auto secondOrderVecVec = secondOrderCentralDifference(vec);

    return vectorVectorToDenseMatrix(secondOrderVecVec);
}

/**
 * Will be removed from the interface soon
 */
double Function::eval(DenseVector x) const
{
    auto vec = denseVectorToVector(x);

    return eval(vec);
}

DenseMatrix Function::evalJacobian(DenseVector x) const
{
    auto vec = denseVectorToVector(x);

    auto jacobian = evalJacobian(vec);

    return vectorToDenseVector(jacobian);
}

DenseMatrix Function::evalHessian(DenseVector x) const
{
    auto vec = denseVectorToVector(x);

    auto hessian = evalHessian(vec);

    return vectorVectorToDenseMatrix(hessian);
}

DenseMatrix Function::centralDifference(DenseVector x) const
{
    DenseMatrix dx(1, x.size());

    double h = 1e-6; // perturbation step size
    double hForward = 0.5*h;
    double hBackward = 0.5*h;

    for (unsigned int i = 0; i < getNumVariables(); ++i)
    {
        DenseVector xForward(x);
        xForward(i) = xForward(i) + hForward;

        DenseVector xBackward(x);
        xBackward(i) = xBackward(i) - hBackward;

        double yForward = eval(xForward);
        double yBackward = eval(xBackward);

        dx(i) = (yForward - yBackward)/(hBackward + hForward);
    }

    return dx;
}

} // namespace SPLINTER
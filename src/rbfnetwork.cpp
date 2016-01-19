/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "rbfnetwork.h"
#include <serializer.h>
#include <utilities.h>
#include "linearsolvers.h"
#include "Eigen/SVD"


namespace SPLINTER
{

RBFNetwork::RBFNetwork(const char *fileName)
    : RBFNetwork(std::string(fileName))
{
}

RBFNetwork::RBFNetwork(const std::string &fileName)
    : LinearFunction(1, DenseVector::Zero(1))
{
    load(fileName);
}

DenseVector RBFNetwork::evalBasis(DenseVector x) const
{
    checkInput(x);

    auto xv = denseVectorToVector(x);

    int i = 0;
    DenseVector basis(getNumCoefficients());
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
        basis(i) = fn->eval(dist(xv, it->getX()));

    if (normalized) basis /= basis.sum();

    return basis;
}

DenseMatrix RBFNetwork::evalJacobian(DenseVector x) const
{
    checkInput(x);

    auto x_vec = denseVectorToVector(x);

    DenseMatrix jac = DenseMatrix::Zero(1, numVariables);

    for (unsigned int i = 0; i < numVariables; i++)
    {
        double sumw = 0;
        double sumw_d = 0;
        double sum = 0;
        double sum_d = 0;

        int j = 0;
        for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++j)
        {
            // Sample
            auto s_vec = it->getX();

            // Distance from sample
            double r = dist(x_vec, s_vec);
            double ri = x_vec.at(i) - s_vec.at(i);

            // Evaluate RBF and its derivative at r
            double f = fn->eval(r);
            double dfdr = fn->evalDerivative(r);

            sum += f;
            sumw += coefficients(j) * f;

            // TODO: check if this assumption is correct
            if (r != 0)
            {
                sum_d += dfdr*ri/r;
                sumw_d += coefficients(j) * dfdr * ri / r;
            }
        }

        if (normalized)
            jac(i) = (sum*sumw_d - sum_d*sumw)/(sum*sum);
        else
            jac(i) = sumw_d;
    }
    return jac;
}

void RBFNetwork::save(const std::string &fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void RBFNetwork::load(const std::string &fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

std::string RBFNetwork::getDescription() const {
    std::string description("RadialBasisFunction of type ");
    switch(type) {
    case RBFType::GAUSSIAN:
        description.append("Gaussian");
        break;
    case RBFType::INVERSE_MULTIQUADRIC:
        description.append("Inverse multiquadric");
        break;
    case RBFType::INVERSE_QUADRIC:
        description.append("Inverse quadric");
        break;
    case RBFType::MULTIQUADRIC:
        description.append("Multiquadric");
        break;
    case RBFType::THIN_PLATE_SPLINE:
        description.append("Thin plate spline");
        break;
    default:
        description.append("Error: Unknown!");
        break;
    }

    return description;
}

} // namespace SPLINTER

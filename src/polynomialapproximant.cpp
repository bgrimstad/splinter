/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <serializer.h>
#include "polynomialapproximant.h"
#include "linearsolvers.h"
#include "mykroneckerproduct.h"
#include "ols.h"

namespace SPLINTER
{

PolynomialApproximant::PolynomialApproximant()
    : Approximant(1),
      poly(Polynomial(1, 1))
{
}

PolynomialApproximant::PolynomialApproximant(const char *fileName)
    : PolynomialApproximant(std::string(fileName))
{
}

PolynomialApproximant::PolynomialApproximant(const std::string fileName)
    : Approximant(1),
      poly(Polynomial(1, 1))
{
    load(fileName);
}

PolynomialApproximant::PolynomialApproximant(const Sample &sample, unsigned int degree)
    : PolynomialApproximant(sample, std::vector<unsigned int>(sample.getNumVariables(), degree))
{
}

PolynomialApproximant::PolynomialApproximant(const Sample &sample, std::vector<unsigned int> degrees)
    : Approximant(sample.getNumVariables()),
      poly(Polynomial(degrees))
{
    if (poly.getNumCoefficients() > sample.size())
        throw Exception("PolynomialApproximant::PolynomialApproximant: Insufficient number of samples!");

    DenseMatrix coefficients = computeCoefficients(poly, sample);

    poly.setCoefficients(coefficients);
}

double PolynomialApproximant::eval(DenseVector x) const
{
    return poly.eval(x);
}

DenseMatrix PolynomialApproximant::evalJacobian(DenseVector x) const
{
    return poly.evalJacobian(x);
}

void PolynomialApproximant::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void PolynomialApproximant::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string PolynomialApproximant::getDescription() const
{
    std::string description("PolynomialApproximant based on polynomial:");
    description.append(poly.getDescription());
    return description;
}

} // namespace SPLINTER

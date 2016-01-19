/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "polynomialbuilder.h"
#include "ols.h"

namespace SPLINTER
{

Polynomial Polynomial::Builder::build() const
{
    Polynomial poly(_degrees, DenseVector::Zero(computeNumBasisFunctions(_degrees)));

    if (poly.getNumCoefficients() > _data.getNumSamples())
        throw Exception("Polynomial::Builder::build: Insufficient number of samples!");

    poly.numVariables = _data.getNumVariables();
    poly.coefficients = computeCoefficients(poly, _data);
    poly.degrees = _degrees;

    return poly;
}

unsigned int Polynomial::Builder::computeNumBasisFunctions(std::vector<unsigned int> degrees) const
{
    unsigned int numMonomials = 1;
    for (auto deg : degrees)
        numMonomials *= (deg+1);

    return numMonomials;
}

} // namespace SPLINTER

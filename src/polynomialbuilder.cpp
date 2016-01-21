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
    Polynomial poly(_powers);

    if (poly.getNumCoefficients() > _data.getNumSamples())
        throw Exception("Polynomial::Builder::build: Insufficient number of samples!");

    poly.numVariables = _data.getNumVariables();
    poly.coefficients = computeCoefficients(poly, _data);
    poly.powers = _powers;

    return poly;
}

} // namespace SPLINTER

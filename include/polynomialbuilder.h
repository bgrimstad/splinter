/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIALBUILDER_H
#define SPLINTER_POLYNOMIALBUILDER_H

#include "builderbase.h"
#include "polynomial.h"
#include "datatable.h"
#include <iostream>

namespace SPLINTER
{

class SPLINTER_API Polynomial::Builder : public BuilderBase_CRTP<Polynomial>
{
public:
    Builder(const DataTable &data) :
            BuilderBase_CRTP(data),
            _degrees(std::vector<unsigned int>(_data.getNumVariables(), 0))
    {}

    // Set build options
    Builder& degree(std::vector<unsigned int> degrees)
    {
        _degrees = degrees;
        return *this;
    }

    Builder& degree(unsigned int degree)
    {
        _degrees = std::vector<unsigned int>(_data.getNumVariables(), degree);
        return *this;
    }

    // Build RBFNetwork
    Polynomial build() const override;

private:
    Builder();

    unsigned int computeNumBasisFunctions(std::vector<unsigned int> degrees) const;

    std::vector<unsigned int> _degrees;
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIALBUILDER_H

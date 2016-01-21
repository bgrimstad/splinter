/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIAL2BUILDER_H
#define SPLINTER_POLYNOMIAL2BUILDER_H

#include "builderbase.h"
#include "polynomial.h"
#include "datatable.h"

namespace SPLINTER
{

class SPLINTER_API Polynomial::Builder : public BuilderBase_CRTP<Polynomial>
{
public:
    Builder(const DataTable &table)
            :
            BuilderBase_CRTP(table),
            _powers(DenseMatrix::Zero(0, 0))
    {}

    Builder& lambda(double lambda)
    {
        if (lambda < 0)
            throw Exception("Polynomial::Builder::lambda: Lambda must be non-negative.");

        _lambda = lambda;
        return *this;
    }

    // Set build options
    Builder& powers(DenseMatrix powers)
    {
        _powers = powers;
        return *this;
    }

    // Build Polynomial
    Polynomial build() const override;

private:
    Builder();

    DenseMatrix _powers;
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIAL2BUILDER_H

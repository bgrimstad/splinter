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
#include "polynomial2.h"
#include "datatable.h"

namespace SPLINTER
{

class SPLINTER_API Polynomial2::Builder : public BuilderBase_CRTP<Polynomial2>
{
public:
    Builder(const DataTable &table)
            :
            BuilderBase_CRTP(table),
            _powers(DenseMatrix::Zero(0, 0))
    {}

    // Set build options
    Builder& powers(DenseMatrix powers)
    {
        _powers = powers;
        return *this;
    }

    // Build Polynomial2
    Polynomial2 build() const override;

private:
    Builder();

    DenseMatrix _powers;
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIAL2BUILDER_H

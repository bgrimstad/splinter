/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BUILDERBASE_H
#define SPLINTER_BUILDERBASE_H

#include "datatable.h"
#include <iostream>

namespace SPLINTER
{

/**
 * Pure abstract builder base class.
 * Used for building Functions
 */
template <class T>
class BuilderBase
{
public:
    BuilderBase(const DataTable &data)
        :
        _data(data),
        _lambda(0.0)
    {}
    virtual ~BuilderBase() {};

    virtual T build() const = 0;

    BuilderBase<T> &lambda(double lambda)
    {
        if (lambda < 0)
        {
            throw Exception("Builder::lambda: Lambda must be non-negative.");
        }

        _lambda = lambda;
        return *this;
    }

protected:
    DataTable _data;
    double _lambda;

private:
    BuilderBase();
};

} // namespace SPLINTER

#endif // SPLINTER_BUILDERBASE_H

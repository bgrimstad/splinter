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

#include "function.h"
#include "datatable.h"
#include <iostream>

namespace SPLINTER
{

/**
 * Pure abstract builder base class.
 * Used for building Functions.
 * Using the Curiously recurring template pattern.
 */
template <class T>
class BuilderBase_CRTP
{
public:
    BuilderBase_CRTP(const DataTable &data)
        : _data(data),
         _lambda(0.0)
    {}
    virtual ~BuilderBase_CRTP() {};

    /*
     * Used by the C interface.
     * Uses polymorphism so we need to be able to get a pointer.
     */
    virtual Function *build_ptr() const
    {
        return new T(build());
    }

    virtual T build() const = 0;

    int getNumVariables() const
    {
        return _data.getNumVariables();
    }

protected:
    DataTable _data;
    double _lambda;

private:
    BuilderBase_CRTP();
};

} // namespace SPLINTER

#endif // SPLINTER_BUILDERBASE_H

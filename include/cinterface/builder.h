/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BUILDER_H
#define SPLINTER_BUILDER_H

#include "cinterface/cinterface.h"
#include "cinterface/utilities.h"

namespace SPLINTER
{

template <class T>
void builder_set_lambda(splinter_obj_ptr builder_ptr, double lambda)
{
    typename T::Builder *builder = get_builder<T>(builder_ptr);
    if (builder == nullptr)
    {
        return;
    }

    builder->lambda(lambda);
}

template <class T>
T *builder_build(splinter_obj_ptr builder_ptr)
{
    typename T::Builder *builder = get_builder<T>(builder_ptr);
    if (builder == nullptr)
    {
        return nullptr;
    }

    return new T(builder->build());
}

template <class T>
void builder_delete(splinter_obj_ptr builder_ptr)
{
    typename T::Builder *builder = get_builder<T>(builder_ptr);

    // No need to check for nullptr, deleting a nullptr has no effect
    delete builder;
}

} // namespace SPLINTER

#endif // SPLINTER_BUILDER_H

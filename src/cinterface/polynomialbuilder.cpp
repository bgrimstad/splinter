/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <cinterface/builder.h>
#include "polynomial.h"
#include "polynomialbuilder.h"
#include "cinterface/cinterface.h"
#include "cinterface/utilities.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_polynomial_builder_init(splinter_obj_ptr datatable_ptr)
{
    splinter_obj_ptr polynomial_builder_ptr = nullptr;

    try
    {
        DataTable *dataTable = get_datatable(datatable_ptr);
        polynomial_builder_ptr = new Polynomial::Builder(*dataTable);
        builders.insert(polynomial_builder_ptr);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return polynomial_builder_ptr;
}

void splinter_polynomial_builder_set_degree(splinter_obj_ptr polynomial_builder_ptr, unsigned int *degrees, int n)
{
    auto builder = get_builder<Polynomial>(polynomial_builder_ptr);
    if(builder != nullptr)
    {
        auto _degrees = get_vector(degrees, n);
        builder->degree(_degrees);
    }
}

void splinter_polynomial_builder_set_lambda(splinter_obj_ptr polynomial_builder_ptr, double lambda)
{
    builder_set_lambda<Polynomial>(polynomial_builder_ptr, lambda);
}

splinter_obj_ptr splinter_polynomial_builder_build(splinter_obj_ptr polynomial_builder_ptr)
{
    auto func = builder_build<Polynomial>(polynomial_builder_ptr);

    if (func != nullptr)
    {
        functions.insert(func);
    }

    return func;
}

void splinter_polynomial_builder_delete(splinter_obj_ptr polynomial_builder_ptr)
{
    builder_delete<Polynomial>(polynomial_builder_ptr);
}

} // extern "C"

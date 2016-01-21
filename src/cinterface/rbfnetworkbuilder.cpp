/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <cinterface/builder.h>
#include "rbfnetwork.h"
#include "rbfnetworkbuilder.h"
#include "cinterface/cinterface.h"
#include "cinterface/utilities.h"

using namespace SPLINTER;

extern "C"
{

void splinter_rbfnetwork_builder_set_type(splinter_obj_ptr rbfnetwork_builder_ptr, int type_index)
{
    RBFType type;
    switch (type_index)
    {
        case 1:
            type = RBFType::THIN_PLATE_SPLINE;
            break;
        case 2:
            type = RBFType::MULTIQUADRIC;
            break;
        case 3:
            type = RBFType::INVERSE_QUADRIC;
            break;
        case 4:
            type = RBFType::INVERSE_MULTIQUADRIC;
            break;
        case 5:
            type = RBFType::GAUSSIAN;
            break;
        default:
            type = RBFType::THIN_PLATE_SPLINE;
            break;
    }

    auto builder = get_builder<RBFNetwork>(rbfnetwork_builder_ptr);
    if(builder != nullptr)
    {
        builder->type(type);
    }
}

void splinter_rbfnetwork_builder_set_normalized(splinter_obj_ptr rbfnetwork_builder_ptr, int normalized)
{
    auto builder = get_builder<RBFNetwork>(rbfnetwork_builder_ptr);
    if(builder != nullptr)
    {
        builder->normalized((bool) normalized);
    }
}

void splinter_rbfnetwork_builder_set_precondition(splinter_obj_ptr rbfnetwork_builder_ptr, int precondition)
{
    auto builder = get_builder<RBFNetwork>(rbfnetwork_builder_ptr);
    if(builder != nullptr)
    {
        builder->precondition((bool) precondition);
    }
}

void splinter_rbfnetwork_builder_set_lambda(splinter_obj_ptr rbfnetwork_builder_ptr, double lambda)
{
    builder_set_lambda<RBFNetwork>(rbfnetwork_builder_ptr, lambda);
}

splinter_obj_ptr splinter_rbfnetwork_builder_build(splinter_obj_ptr rbfnetwork_builder_ptr)
{
    auto func = builder_build<RBFNetwork>(rbfnetwork_builder_ptr);

    if (func != nullptr)
    {
        functions.insert(func);
    }

    return func;
}

void splinter_rbfnetwork_builder_delete(splinter_obj_ptr rbfnetwork_builder_ptr)
{
    builder_delete<RBFNetwork>(rbfnetwork_builder_ptr);
}

} // extern "C"

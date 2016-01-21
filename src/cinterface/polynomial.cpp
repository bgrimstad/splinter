/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "cinterface/cinterface.h"
#include "cinterface/utilities.h"
#include "polynomial.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_polynomial_load_init(const char *filename)
{
    splinter_obj_ptr polyfit = nullptr;

    try
    {
        polyfit = (splinter_obj_ptr) new Polynomial(filename);
        functions.insert(polyfit);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return polyfit;
}

} // extern "C"

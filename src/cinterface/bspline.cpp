/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline.h"
#include "cinterface/utilities.h"
#include "cinterface/cinterface.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_bspline_load_init(const char *filename)
{
    splinter_obj_ptr bspline = nullptr;

    try
    {
        bspline = (splinter_obj_ptr) new BSpline(filename);
        functions.insert(bspline);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

} // extern "C"

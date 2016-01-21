/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <cinterface/utilities.h>
#include "rbfnetwork.h"
#include "cinterface/cinterface.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_rbf_load_init(const char *filename)
{
    splinter_obj_ptr rbf = nullptr;

    try
    {
        rbf = (splinter_obj_ptr) new RBFNetwork(filename);
        functions.insert(rbf);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return rbf;
}

} // extern "C"

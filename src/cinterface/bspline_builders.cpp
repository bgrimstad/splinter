/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_builders.h"
#include "cinterface/cinterface.h"
#include "cinterface/utilities.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_bspline_interpolator(splinter_obj_ptr datatable_ptr, int degree)
{
    splinter_obj_ptr bspline = nullptr;
    DataTable *data_table = get_datatable(datatable_ptr);
    try
    {
        bspline = (splinter_obj_ptr) bspline_interpolator(*data_table, degree).clone();  // TODO: unnecessary copy
        bsplines.insert(bspline);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

splinter_obj_ptr splinter_bspline_smoother(splinter_obj_ptr datatable_ptr, int degree, int smoothing, double alpha,
                                           double *weights, unsigned int num_weights)
{
    splinter_obj_ptr bspline = nullptr;

    try
    {
        DataTable *data_table = get_datatable(datatable_ptr);
        auto _degree = static_cast<unsigned int>(degree);
        auto _weights = get_vector(weights, num_weights);
        auto _smoothing = resolve_smoothing(smoothing);

        // TODO: Fix unnecessary copy
        bspline = (splinter_obj_ptr) bspline_smoother(*data_table, _degree, _smoothing, alpha, _weights).clone();
        bsplines.insert(bspline);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

splinter_obj_ptr splinter_bspline_unfitted(splinter_obj_ptr datatable_ptr, unsigned int *degrees,
                                           unsigned int num_degrees, int knot_spacing,
                                           unsigned int *num_basis_functions,
                                           unsigned int num_num_basis_functions)
{
    splinter_obj_ptr bspline = nullptr;

    try
    {
        DataTable *data_table = get_datatable(datatable_ptr);
        auto _degrees = get_vector(degrees, num_degrees);
        auto _knot_spacing = resolve_knot_spacing(knot_spacing);
        auto _num_basis_functions = get_vector(num_basis_functions, num_num_basis_functions);

        // TODO: Fix unnecessary copy
        bspline = (splinter_obj_ptr) bspline_unfitted(*data_table, _degrees, _knot_spacing, _num_basis_functions).clone();
        bsplines.insert(bspline);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;

}

} // extern "C"

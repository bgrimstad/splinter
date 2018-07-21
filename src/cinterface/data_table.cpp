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
#include "data_table.h"

// Used for printing debug info to files when the library is loaded from eg. Python
#ifndef NDEBUG
# include <fstream>
#endif

using namespace SPLINTER;

extern "C"
{

/* DataTable constructor */
splinter_obj_ptr splinter_datatable_init()
{
    splinter_obj_ptr dataTable = (splinter_obj_ptr) new DataTable();

    datatables.insert(dataTable);

    return dataTable;
}

void splinter_datatable_add_samples_row_major(splinter_obj_ptr datatable_ptr,
                                              double *xs, int x_dim,
                                              double *ys, int y_dim,
                                              int n_samples)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        try
        {
            std::vector<double> x_vec(x_dim, 0);
            std::vector<double> y_vec(y_dim, 0);

            for (int i = 0; i < n_samples; ++i)
            {
                // Vectors have been initialised to appropriate sizes, so this should be safe
                memcpy(x_vec.data(), &xs[x_dim*i], sizeof(double) * x_dim);
                memcpy(y_vec.data(), &ys[y_dim*i], sizeof(double) * y_dim);

                dataTable->add_sample(x_vec, y_vec);
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void splinter_datatable_add_samples_col_major(splinter_obj_ptr datatable_ptr, double *x, int n_samples, int x_dim)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        try
        {
            std::vector<double> vec(x_dim, 0);
            for (int i = 0; i < n_samples; ++i)
            {
                for (int j = 0; j < x_dim; ++j)
                {
                    vec.at(j) = x[i + j * n_samples];
                }

                dataTable->add_sample(vec, x[i + x_dim * n_samples]);
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

int splinter_datatable_get_dim_x(splinter_obj_ptr datatable_ptr)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return (int) dataTable->get_dim_x();
    }

    return 0;
}

int splinter_datatable_get_dim_y(splinter_obj_ptr datatable_ptr)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return (int) dataTable->get_dim_y();
    }

    return 0;
}

int splinter_datatable_get_num_samples(splinter_obj_ptr datatable_ptr)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return dataTable->get_num_samples();
    }

    return 0;
}

void splinter_datatable_to_json(splinter_obj_ptr datatable_ptr, const char *filename)
{
    auto datatable = get_datatable(datatable_ptr);
    if (datatable != nullptr)
    {
        try
        {
            datatable->to_json(filename);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

splinter_obj_ptr splinter_datatable_from_json(const char *filename)
{
    splinter_obj_ptr loaded_datatable = nullptr;

    try
    {
        loaded_datatable = (splinter_obj_ptr) new DataTable(DataTable::from_json(filename));
        datatables.insert(loaded_datatable);
    }
    catch (const Exception &e)
    {
        // Free the memory of the copy in case the exception was thrown by the insert
        delete (DataTable *) loaded_datatable;
        set_error_string(e.what());
    }

    return loaded_datatable;
}

void splinter_datatable_delete(splinter_obj_ptr datatable_ptr)
{
    auto dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        datatables.erase(datatable_ptr);
        delete dataTable;
    }
}

} // extern "C"

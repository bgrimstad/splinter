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

using namespace SPLINTER;

extern "C"
{

double *splinter_function_eval_row_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        try
        {
            size_t num_variables = func->getNumVariables();
            size_t num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector<double>(x, num_variables);
                retVal[i] = func->eval(xvec);
                x += num_variables;
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return retVal;
}

double *splinter_function_eval_jacobian_row_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        try
        {
            size_t num_variables = func->getNumVariables();
            size_t num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_variables * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector<double>(x, num_variables);
                DenseMatrix jacobian = func->evalJacobian(xvec);

                /* Copy jacobian from stack to heap */
                memcpy(retVal + i*num_variables, jacobian.data(), sizeof(double) * num_variables);
                x += num_variables;
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return retVal;
}

double *splinter_function_eval_hessian_row_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        try
        {
            size_t num_variables = func->getNumVariables();
            size_t num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_variables * num_variables * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector<double>(x, num_variables);
                DenseMatrix hessian = func->evalHessian(xvec);

                /* Copy hessian from stack to heap */
                memcpy(retVal + i*num_variables*num_variables, hessian.data(), sizeof(double) * num_variables * num_variables);
                x += num_variables;
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return retVal;
}

double *splinter_function_eval_col_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, func->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = splinter_function_eval_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *splinter_function_eval_jacobian_col_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, func->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = splinter_function_eval_jacobian_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *splinter_function_eval_hessian_col_major(splinter_obj_ptr function, double *x, int x_len)
{
    double *retVal = nullptr;

    auto func = get_function(function);
    if (func != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, func->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = splinter_function_eval_hessian_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }
    return retVal;
}

int splinter_function_get_num_variables(splinter_obj_ptr function)
{
    int retVal = 0;

    auto func = get_function(function);
    if (func != nullptr)
    {
        retVal = func->getNumVariables();
    }

    return retVal;
}

void splinter_function_save(splinter_obj_ptr function, const char *filename)
{
    auto func = get_function(function);
    if (func != nullptr)
    {
        try
        {
            func->save(filename);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void splinter_function_delete(splinter_obj_ptr function)
{
    auto func = get_function(function);

    if (func != nullptr)
    {
        functions.erase(function);
        delete func;
    }
}

} // extern "C"

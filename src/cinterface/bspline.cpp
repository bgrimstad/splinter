/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <fstream>
#include "bspline.h"
#include "cinterface/utilities.h"

using namespace SPLINTER;

extern "C"
{

splinter_obj_ptr splinter_bspline_from_param(unsigned int dim_x, unsigned int dim_y, unsigned int *degrees,
                                             double *knot_vectors, unsigned int *num_knots_per_vector,
                                             double *control_points, unsigned int num_control_points)
{
    splinter_obj_ptr bspline = nullptr;
    auto degrees_vec = get_vector(degrees, dim_x);
    auto knot_vectors_vec_vec = get_vector_vector(knot_vectors, num_knots_per_vector, dim_x);
    auto control_points_vec_vec = get_vector_vector(control_points, dim_y, num_control_points);

    try
    {
        bspline = (splinter_obj_ptr) new BSpline(degrees_vec, knot_vectors_vec_vec, control_points_vec_vec);
        bsplines.insert(bspline);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

splinter_obj_ptr splinter_bspline_from_param_zero(unsigned int dim_x,
                                                  unsigned int dim_y,
                                                  unsigned int *degrees,
                                                  double *knot_vectors,
                                                  unsigned int *num_knots_per_vector)
{
    splinter_obj_ptr bspline = nullptr;
    auto knot_vectors_vec_vec = get_vector_vector(knot_vectors, num_knots_per_vector, dim_x);
    auto degrees_vec = get_vector(degrees, dim_x);

    try
    {
        bspline = (splinter_obj_ptr) new BSpline(degrees_vec, knot_vectors_vec_vec, dim_y);
        bsplines.insert(bspline);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

int *splinter_bspline_get_knot_vector_sizes(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    int *sizes = nullptr;
    if (bspline != nullptr)
    {
        try
        {
            auto knot_vectors = bspline->get_knot_vectors();

            sizes = (int *) malloc(knot_vectors.size() * sizeof (int));

            if (sizes != nullptr)
            {
                int i = 0;
                for (auto knot_vector : knot_vectors)
                {
                    sizes[i++] = (int) knot_vector.size();
                }
            }
            else
            {
                set_error_string("Unable to allocate memory!");
            }
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return sizes;
}

double *splinter_bspline_get_knot_vectors(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    double *knot_vectors_as_array = nullptr;
    if (bspline != nullptr)
    {
        try
        {
            auto knot_vectors = bspline->get_knot_vectors();

            // Overkill, but some C++11 is nice
            int total_n_elements = 0;
            std::for_each(knot_vectors.cbegin(), knot_vectors.cend(),
                          [&total_n_elements](const std::vector<double> &knot_vector)
                          {
                              total_n_elements += knot_vector.size();
                          });
            knot_vectors_as_array = (double *) malloc(total_n_elements * sizeof (double));

            if (knot_vectors_as_array != nullptr)
            {
                int i = 0;
                for (auto knot_vector : knot_vectors)
                {
                    // Even more unnecessary C++11 stuff
                    std::copy(knot_vector.begin(), knot_vector.end(), knot_vectors_as_array + i);
                    i += knot_vector.size();
                }
            }
            else
            {
                set_error_string("Unable to allocate memory!");
            }
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return knot_vectors_as_array;
}

int splinter_bspline_get_num_control_points(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            return bspline->get_num_control_points();
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return -1;
}

double *splinter_bspline_get_control_points(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    double *control_points_as_array = nullptr;
    if (bspline != nullptr)
    {
        try
        {
            auto control_points = bspline->get_control_points();

            control_points_as_array = (double *) malloc(control_points.size() * sizeof (double));

            if (control_points_as_array != nullptr)
            {
                for (int i = 0; i < control_points.rows(); ++i)
                {
                    for (int j = 0; j < control_points.cols(); ++j)
                    {
                        control_points_as_array[i*control_points.cols() + j] = control_points(i, j);
                    }
                }
            }
            else
            {
                set_error_string("Unable to allocate memory!");
            }
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return control_points_as_array;
}

double *splinter_bspline_get_knot_averages(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    double *knot_averages_as_array = nullptr;
    if (bspline != nullptr)
    {
        try
        {
            auto knot_averages = bspline->get_knot_averages();

            knot_averages_as_array = (double *) malloc(knot_averages.size() * sizeof (double));

            if (knot_averages_as_array != nullptr)
            {
                for (int i = 0; i < knot_averages.rows(); ++i)
                {
                    for (int j = 0; j < knot_averages.cols(); ++j)
                    {
                        knot_averages_as_array[i*knot_averages.cols() + j] = knot_averages(i, j);
                    }
                }
            }
            else
            {
                set_error_string("Unable to allocate memory!");
            }
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return knot_averages_as_array;
}

int *splinter_bspline_get_basis_degrees(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    int *basis_degrees_as_array = nullptr;
    if (bspline != nullptr)
    {
        try
        {
            auto basis_degrees = bspline->get_basis_degrees();

            basis_degrees_as_array = (int *) malloc(basis_degrees.size() * sizeof (int));

            if (basis_degrees_as_array != nullptr)
            {
                for (unsigned int i = 0; i < basis_degrees.size(); ++i)
                {
                    basis_degrees_as_array[i] = basis_degrees[i];
                }
            }
            else
            {
                set_error_string("Unable to allocate memory!");
            }
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
    return basis_degrees_as_array;
}

double *splinter_bspline_eval_row_major(splinter_obj_ptr bspline_ptr, double *xs, int x_len)
{
    double *retVal = nullptr;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            size_t x_dim = bspline->get_dim_x();
            size_t y_dim = bspline->get_dim_y();
            size_t num_points = x_len / x_dim;

            retVal = (double *) malloc(sizeof(double) * num_points * y_dim);
            for (size_t i = 0; i < num_points; i++)
            {
                auto xvec = get_vector<double>(xs, x_dim);
                // Underlying memory of vector is guaranteed to be contiguous
                memcpy(&retVal[i*y_dim], bspline->eval(xvec).data(), sizeof(double) * y_dim);
                xs += x_dim;
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return retVal;
}

double *splinter_bspline_eval_jacobian_row_major(splinter_obj_ptr bspline_ptr, double *x, int x_len)
{
    double *retVal = nullptr;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            size_t num_variables = bspline->get_dim_x();
            size_t num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_variables * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector<double>(x, num_variables);
                DenseMatrix jacobian = bspline->eval_jacobian(xvec);

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

double *splinter_bspline_eval_col_major(splinter_obj_ptr bspline_ptr, double *x, int x_len)
{
    double *retVal = nullptr;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, bspline->get_dim_x(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = splinter_bspline_eval_row_major(bspline, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *splinter_bspline_eval_jacobian_col_major(splinter_obj_ptr bspline_ptr, double *x, int x_len)
{
    double *retVal = nullptr;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, bspline->get_dim_x(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = splinter_bspline_eval_jacobian_row_major(bspline, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

int splinter_bspline_get_dim_x(splinter_obj_ptr bspline_ptr)
{
    int retVal = 1;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        retVal = bspline->get_dim_x();
    }

    return retVal;
}

int splinter_bspline_get_dim_y(splinter_obj_ptr bspline_ptr)
{
    int retVal = 1;

    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        retVal = bspline->get_dim_y();
    }

    return retVal;
}

void splinter_bspline_to_json(splinter_obj_ptr bspline_ptr, const char *filename)
{
    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            bspline->to_json(filename);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

splinter_obj_ptr splinter_bspline_from_json(const char *filename)
{
    splinter_obj_ptr loaded_bspline = nullptr;

    try
    {
        loaded_bspline = (splinter_obj_ptr) new BSpline(BSpline::from_json(filename));
        bsplines.insert(loaded_bspline);
    }
    catch (const Exception &e)
    {
        // Free the memory of the copy in case the exception was thrown by the insert
        delete (BSpline *) loaded_bspline;
        set_error_string(e.what());
    }

    return loaded_bspline;
}

void splinter_bspline_delete(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);

    if (bspline != nullptr)
    {
        bsplines.erase(bspline_ptr);
        delete bspline;
    }
}

void splinter_bspline_insert_knots(splinter_obj_ptr bspline_ptr, double tau, unsigned int dim, unsigned int multiplicity)
{
    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            bspline->insert_knots(tau, dim, multiplicity);
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void splinter_bspline_decompose_to_bezier_form(splinter_obj_ptr bspline_ptr)
{
    auto bspline = get_bspline(bspline_ptr);
    if (bspline != nullptr)
    {
        try
        {
            bspline->decompose_to_bezier();
        }
        catch (const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

splinter_obj_ptr splinter_bspline_copy(splinter_obj_ptr bspline_ptr) {
    auto bspline = get_bspline(bspline_ptr);
    splinter_obj_ptr copy = nullptr;

    if (bspline != nullptr)
    {
        try
        {
            copy = (splinter_obj_ptr) bspline->clone();
            bsplines.insert(copy);
        }
        catch (const Exception &e)
        {
            // Free the memory of the copy in case the exception was thrown by the insert
            delete (BSpline *) copy;
            set_error_string(e.what());
        }
    }
    return copy;
}

splinter_obj_ptr splinter_bspline_fit(splinter_obj_ptr bspline_ptr, splinter_obj_ptr datatable_ptr, int smoothing,
                                      double alpha, double *weights, int num_weights)
{
    auto bspline = get_bspline(bspline_ptr);
    if (bspline == nullptr)
    {
        return nullptr;
    }

    try
    {
        DataTable *dataTable = get_datatable(datatable_ptr);
        auto _smoothing = resolve_smoothing(smoothing);
        auto _weights = get_vector(weights, num_weights);
        auto new_bspline = bspline->clone();
        new_bspline->fit(*dataTable, _smoothing, alpha, _weights);
        bsplines.insert(new_bspline);
        return new_bspline;
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return nullptr;
}

} // extern "C"

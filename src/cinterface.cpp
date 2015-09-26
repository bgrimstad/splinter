/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "cinterface.h"
#include <sample.h>
#include <bspline.h>
#include <bsplineapproximant.h>
#include <psplineapproximant.h>
#include <rbfapproximant.h>
#include <polynomialapproximant.h>
#include "definitions.h"
#include <set>
#include <iostream>

using namespace SPLINTER;

// 1 if the last function call caused an error, 0 else
int lastFuncCallError = 0;

const char *error_string = "No error.";

// Keep a list of objects so we avoid performing operations on objects that don't exist
std::set<obj_ptr> objects = std::set<obj_ptr>();

static void set_error_string(const char *new_error_string)
{
    error_string = new_error_string;
    lastFuncCallError = 1;
}

/* Cast the obj_ptr to a DataTable * */
static Sample *get_datatable(obj_ptr datatable_ptr)
{
    if (objects.count(datatable_ptr) > 0)
    {
        return (Sample *) datatable_ptr;
    }

    set_error_string("Invalid reference to DataTable: Maybe it has been deleted?");

    return nullptr;
}

/* Cast the obj_ptr to an Approximant * */
static Approximant *get_approximant(obj_ptr approximant_ptr)
{
    if (objects.count(approximant_ptr) > 0)
    {
        return (Approximant *) approximant_ptr;
    }

    set_error_string("Invalid reference to Approximant: Maybe it has been deleted?");

    return nullptr;
}

static DenseVector get_densevector(double *x, int x_dim)
{
    DenseVector xvec(x_dim);
    for (int i = 0; i < x_dim; i++)
    {
        xvec(i) = x[i];
    }

    return xvec;
}

static double *get_row_major(double *col_major, size_t point_dim, size_t x_len)
{
    if (point_dim == 0)
    {
        set_error_string("Dimension of x should be larger than 0!");
        return nullptr;
    }

    double *row_major = (double *) malloc(sizeof(double) * x_len);
    if(row_major == nullptr)
    {
        set_error_string("Out of memory!");
        return nullptr;
    }

    size_t num_points = x_len / point_dim;
    for (size_t i = 0; i < x_len; ++i)
    {
        size_t sample_number = i / point_dim; // Intentional integer division
        size_t dimension_number = i % point_dim;
        row_major[i] = col_major[dimension_number * num_points + sample_number];
    }

    return row_major;
}

extern "C"
{
/* 1 if the last call to the library resulted in an error,
    *  0 otherwise. This is used to avoid throwing exceptions across library boundaries,
    *  and we expect the caller to manually check the value of this flag.
    *  We always reset it when it is checked, then we won't have to reset it in the beginning of every
    *  end user visible function.
    */
int get_error()
{
    int temp = lastFuncCallError;
    lastFuncCallError = 0;
    return temp;
}

const char *get_error_string()
{
    return error_string;
}

/* DataTable constructor */
obj_ptr datatable_init()
{
    obj_ptr dataTable = (obj_ptr) new Sample();

    objects.insert(dataTable);

    return dataTable;
}

obj_ptr datatable_load_init(const char *filename)
{
    obj_ptr dataTable = nullptr;

    try
    {
        dataTable = (obj_ptr) new Sample(filename);
        objects.insert(dataTable);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return dataTable;
}

void datatable_add_samples_row_major(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        try
        {
            DenseVector vec(x_dim);
            for (int i = 0; i < n_samples; ++i)
            {
                int sample_start = i*(x_dim+1);
                for (int offset = 0; offset < x_dim; ++offset)
                {
                    vec(offset) = x[sample_start + offset];
                }

//            std::cout << "Adding sample: (";
//            for(int l = 0; l < vec.size(); l++)
//            {
//                if(l != 0)
//                    std::cout << ",";
//                std::cout << vec(l);
//            }
//            std::cout << ") = ";
//            std::cout << x[sample_start + x_dim] << std::endl;

                dataTable->addSample(vec, x[sample_start + x_dim]);
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void datatable_add_samples_col_major(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim, int size)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        try
        {
            DenseVector vec(x_dim);
            for (int i = 0; i < n_samples; ++i)
            {
                for (int j = 0; j < x_dim; ++j)
                {
                    vec(j) = x[i + j * size];
                }

//            std::cout << "Adding sample: (";
//            for(int l = 0; l < vec.size(); l++)
//            {
//                if(l != 0)
//                    std::cout << ",";
//                std::cout << vec(l);
//            }
//            std::cout << ") = ";
//            std::cout << x[i + x_dim * size] << std::endl;

                dataTable->addSample(vec, x[i + x_dim * size]);
            }
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

int datatable_get_num_variables(obj_ptr datatable_ptr)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return dataTable->getNumVariables();
    }

    return 0;
}


int datatable_get_num_samples(obj_ptr datatable_ptr)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return dataTable->getNumSamples();
    }

    return 0;
}

void datatable_save(obj_ptr datatable_ptr, const char *filename)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        try
        {
            dataTable->save(filename);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void datatable_delete(obj_ptr datatable_ptr)
{
    Sample *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        objects.erase(datatable_ptr);
        delete dataTable;
    }
}

/* BSpline constructor */
obj_ptr bspline_init(obj_ptr datatable_ptr, int degree)
{
    obj_ptr bspline = nullptr;

    auto table = get_datatable(datatable_ptr);
    if (table != nullptr)
    {
        BSplineType bsplineType;
        switch (degree) {
            case 1: {
                bsplineType = BSplineType::LINEAR;
                break;
            }
            case 2: {
                bsplineType = BSplineType::QUADRATIC;
                break;
            }
            case 3: {
                bsplineType = BSplineType::CUBIC;
                break;
            }
            case 4: {
                bsplineType = BSplineType::QUARTIC;
                break;
            }
            default: {
                set_error_string("Invalid BSplineType!");
                return nullptr;
            }
        }

        try
        {
            bspline = (obj_ptr) new BSplineApproximant(*table, bsplineType);
            objects.insert(bspline);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return bspline;
}

obj_ptr bspline_load_init(const char *filename)
{
    obj_ptr bspline = nullptr;

    try
    {
        bspline = (obj_ptr) new BSplineApproximant(filename);
        objects.insert(bspline);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
}

/* PSpline constructor */
obj_ptr pspline_init(obj_ptr datatable_ptr, double lambda)
{
    obj_ptr pspline = nullptr;

    auto table = get_datatable(datatable_ptr);
    if (table != nullptr)
    {
        try
        {
            pspline = (obj_ptr) new PSplineApproximant(*table, lambda);
            objects.insert(pspline);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return pspline;
}

obj_ptr pspline_load_init(const char *filename)
{
    obj_ptr pspline = nullptr;

    try
    {
        pspline = (obj_ptr) new PSplineApproximant(filename);
        objects.insert(pspline);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return pspline;
}

/* RadialBasisFunction constructor */
obj_ptr rbf_init(obj_ptr datatable_ptr, int type_index, int normalized)
{
    obj_ptr rbf = nullptr;

    auto table = get_datatable(datatable_ptr);
    if (table != nullptr)
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

        bool norm = normalized != 0;

        try
        {
            rbf = (obj_ptr) new RBFApproximant(*table, type, norm);
            objects.insert(rbf);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return rbf;
}

obj_ptr rbf_load_init(const char *filename)
{
    obj_ptr rbf = nullptr;

    try
    {
        rbf = (obj_ptr) new RBFApproximant(filename);
        objects.insert(rbf);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return rbf;
}

/* PolynomialApproximant constructor */
obj_ptr polynomial_regression_init(obj_ptr datatable_ptr, int *degrees, int degrees_dim)
{
    obj_ptr polyfit = nullptr;

    auto table = get_datatable(datatable_ptr);
    if (table != nullptr)
    {
        auto degreeVec = std::vector<unsigned int>(degrees_dim);
        for (int i = 0; i < degrees_dim; ++i)
        {
            degreeVec.at(i) = (unsigned int) degrees[i];
        }

        try
        {
            polyfit = (obj_ptr) new PolynomialApproximant(*table, degreeVec);
            objects.insert(polyfit);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }

    return polyfit;
}

obj_ptr polynomial_regression_load_init(const char *filename)
{
    obj_ptr polyfit = nullptr;

    try
    {
        polyfit = (obj_ptr) new PolynomialApproximant(filename);
        objects.insert(polyfit);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return polyfit;
}


double *eval_row_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        try
        {
            int num_variables = approx->getNumVariables();
            int num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector(x, num_variables);
                retVal[i] = approx->eval(xvec);
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

double *eval_jacobian_row_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        try
        {
            int num_variables = approx->getNumVariables();
            int num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_variables * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector(x, num_variables);
                DenseMatrix jacobian = approx->evalJacobian(xvec);

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

double *eval_hessian_row_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        try
        {
            int num_variables = approx->getNumVariables();
            int num_points = x_len / num_variables;

            retVal = (double *) malloc(sizeof(double) * num_variables * num_variables * num_points);
            for (size_t i = 0; i < num_points; ++i)
            {
                auto xvec = get_densevector(x, num_variables);
                DenseMatrix hessian = approx->evalHessian(xvec);

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

double *eval_col_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, approx->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = eval_row_major(approx, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *eval_jacobian_col_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, approx->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = eval_jacobian_row_major(approx, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *eval_hessian_col_major(obj_ptr approximant, double *x, int x_len)
{
    double *retVal = nullptr;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {double *row_major = nullptr;
        try
        {
            row_major = get_row_major(x, approx->getNumVariables(), x_len);
            if (row_major == nullptr)
            {
                return nullptr; // Pass on the error message set by get_row_major
            }

            retVal = eval_hessian_row_major(approx, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }
    return retVal;
}

int approximant_get_num_variables(obj_ptr approximant)
{
    int retVal = 0;

    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        retVal = approx->getNumVariables();
    }

    return retVal;
}

void approximant_save(obj_ptr approximant, const char *filename)
{
    auto approx = get_approximant(approximant);
    if (approx != nullptr)
    {
        try
        {
            approx->save(filename);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
    }
}

void approximant_delete(obj_ptr approximant)
{
    auto approx = get_approximant(approximant);

    if (approx != nullptr)
    {
        objects.erase(approximant);
        delete approx;
    }
}

}

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "cinterface.h"
#include <datatable.h>
#include <bspline.h>
#include <bsplinebuilder.h>
#include "rbfnetwork.h"
#include <rbfbuilder.h>
#include <polynomial.h>
#include <polynomialbuilder.h>
#include "definitions.h"
#include <set>
#include <iostream>

using namespace SPLINTER;

// 1 if the last function call caused an error, 0 else
int lastFuncCallError = 0;

const char *error_string = "No error.";

// Keep a list of objects so we avoid performing operations on objects that don't exist
std::set<obj_ptr> dataTables = std::set<obj_ptr>();
std::set<obj_ptr> functions = std::set<obj_ptr>();
std::set<obj_ptr> bsplineBuilders = std::set<obj_ptr>();

static void set_error_string(const char *new_error_string)
{
    error_string = new_error_string;
    lastFuncCallError = 1;
}

/* Cast the obj_ptr to a DataTable * */
static DataTable *get_datatable(obj_ptr datatable_ptr)
{
    if (dataTables.count(datatable_ptr) > 0)
    {
        return static_cast<DataTable *>(datatable_ptr);
    }

    set_error_string("Invalid reference to DataTable: Maybe it has been deleted?");

    return nullptr;
}

/* Cast the obj_ptr to a Function * */
static Function *get_function(obj_ptr function_ptr)
{
    if (functions.count(function_ptr) > 0)
    {
        return static_cast<Function *>(function_ptr);
    }

    set_error_string("Invalid reference to Function: Maybe it has been deleted?");

    return nullptr;
}

/**
 * Cast bspline_builder_ptr to BSpline::Builder* if possible.
 * Checks the internal objects index to see if bspline_builder_ptr is a valid pointer.
 */
static BSpline::Builder *get_bspline_builder(obj_ptr bspline_builder_ptr)
{
    if (bsplineBuilders.count(bspline_builder_ptr) > 0)
    {
        return static_cast<BSpline::Builder *>(bspline_builder_ptr);
    }

    set_error_string("Invalid reference to BSpline::Builder: Maybe it has been deleted?");

    return nullptr;
}

/**
 * Convert from standard C array to DenseVector.
 *
 * @param x C array to convert from.
 * @param x_dim The size of x.
 * @return DenseVector with the same data as x.
 */
template <class NUMERICAL_TYPE>
static DenseVector get_densevector(NUMERICAL_TYPE *x, size_t x_dim)
{
    DenseVector xvec(x_dim);
    for (size_t i = 0; i < x_dim; i++)
    {
        xvec(i) = (double) x[i];
    }

    return xvec;
}

/**
 * Convert from DenseVector to a vector of NUMERICAL_TYPE.
 * It must be possible to cast from double to NUMERICAL_TYPE.
 *
 * @param x DenseVector to convert from.
 * @return Vector with the same data as x.
 */
template <class NUMERICAL_TYPE>
static std::vector<NUMERICAL_TYPE> get_vector(DenseVector x)
{
    auto vector = std::vector<NUMERICAL_TYPE>(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        vector.at(i) = (NUMERICAL_TYPE) x(i);
    }

    return vector;
}

/**
 * Convert from column major to row major with point_dim number of columns.
 *
 * @param col_major Column major data
 * @param point_dim Dimension of each point (= number of columns)
 * @return col_major data stored row major.
 */
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
        size_t row_num = i / point_dim; // Intentional integer division
        size_t col_num = i % point_dim;
        row_major[i] = col_major[col_num * num_points + row_num];
    }

    return row_major;
}

extern "C"
{
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
    obj_ptr dataTable = (obj_ptr) new DataTable();

    dataTables.insert(dataTable);

    return dataTable;
}

obj_ptr datatable_load_init(const char *filename)
{
    obj_ptr dataTable = nullptr;

    try
    {
        dataTable = (obj_ptr) new DataTable(filename);
        dataTables.insert(dataTable);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return dataTable;
}

void datatable_add_samples_row_major(obj_ptr datatable_ptr, double *x, int n_samples, int x_dim)
{
    DataTable *dataTable = get_datatable(datatable_ptr);
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
    DataTable *dataTable = get_datatable(datatable_ptr);
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
    DataTable *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return dataTable->getNumVariables();
    }

    return 0;
}


int datatable_get_num_samples(obj_ptr datatable_ptr)
{
    DataTable *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        return dataTable->getNumSamples();
    }

    return 0;
}

void datatable_save(obj_ptr datatable_ptr, const char *filename)
{
    DataTable *dataTable = get_datatable(datatable_ptr);
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
    DataTable *dataTable = get_datatable(datatable_ptr);
    if (dataTable != nullptr)
    {
        dataTables.erase(datatable_ptr);
        delete dataTable;
    }
}

obj_ptr bspline_load_init(const char *filename)
{
    obj_ptr bspline = nullptr;

    try
    {
        bspline = (obj_ptr) new BSpline(filename);
        functions.insert(bspline);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline;
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
            RBFNetwork temp = RBFNetwork::Builder(*table).type(type).normalized(norm).build();
            rbf = (obj_ptr) new RBFNetwork(temp);
//            rbf = (obj_ptr) new RBFNetwork(*table, type, norm);
            functions.insert(rbf);
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
        rbf = (obj_ptr) new RBFNetwork(filename);
        functions.insert(rbf);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return rbf;
}

/* Polynomial constructor */
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
            Polynomial temp = Polynomial::Builder(*table).degree(degreeVec).build();
            polyfit = (obj_ptr) new Polynomial(temp);
//            polyfit = (obj_ptr) new Polynomial(*table, degreeVec);
            functions.insert(polyfit);
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
        polyfit = (obj_ptr) new Polynomial(filename);
        functions.insert(polyfit);
    }
    catch(const Exception &e)
    {
        set_error_string(e.what());
    }

    return polyfit;
}

obj_ptr bspline_builder_init(obj_ptr datatable_ptr)
{
    obj_ptr bspline_builder_ptr = nullptr;

    try
    {
        DataTable *dataTable = get_datatable(datatable_ptr);
        bspline_builder_ptr = new BSpline::Builder(*dataTable);
        bsplineBuilders.insert(bspline_builder_ptr);
    }
    catch (const Exception &e)
    {
        set_error_string(e.what());
    }

    return bspline_builder_ptr;
}

void bspline_builder_set_degree(obj_ptr bspline_builder_ptr, int *degrees, int n)
{
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr)
    {
        std::vector<unsigned int> _degrees((unsigned int) n);
        for (int i = 0; i < n; ++i)
        {
            _degrees.at(i) = (unsigned int) degrees[i];
        }
        builder->degree(_degrees);
    }
}

void bspline_builder_set_num_basis_functions(obj_ptr bspline_builder_ptr, int *num_basis_functions, int n)
{
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr)
    {
        std::vector<unsigned int> _num_basis_functions((unsigned int) n);
        for (int i = 0; i < n; ++i)
        {
            _num_basis_functions.at(i) = (unsigned int) num_basis_functions[i];
        }
        builder->numBasisFunctions(_num_basis_functions);
    }
}

void bspline_builder_set_knot_spacing(obj_ptr bspline_builder_ptr, int knot_spacing)
{
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr)
    {
        switch (knot_spacing)
        {
            case 0:
                builder->knotSpacing(BSpline::KnotSpacing::SAMPLE);
                break;
            case 1:
                builder->knotSpacing(BSpline::KnotSpacing::EQUIDISTANT);
                break;
            case 2:
                builder->knotSpacing(BSpline::KnotSpacing::EXPERIMENTAL);
                break;
            default:
                set_error_string("Error: Invalid knot spacing!");
                break;
        }
    }
}

void bspline_builder_set_smoothing(obj_ptr bspline_builder_ptr, int smoothing)
{
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr)
    {
        switch (smoothing)
        {
            case 0:
                builder->smoothing(BSpline::Smoothing::NONE);
                break;
            case 1:
                builder->smoothing(BSpline::Smoothing::REGULARIZATION);
                break;
            case 2:
                builder->smoothing(BSpline::Smoothing::PSPLINE);
                break;
            default:
                set_error_string("Error: Invalid smoothing!");
                break;
        }
    }
}

void bspline_builder_set_lambda(obj_ptr bspline_builder_ptr, double lambda)
{
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr) {
        builder->lambda(lambda);
    }
}

obj_ptr bspline_builder_build(obj_ptr bspline_builder_ptr)
{
    BSpline *bspline = nullptr;
    BSpline::Builder *builder = get_bspline_builder(bspline_builder_ptr);
    if(builder != nullptr)
    {
        bspline = builder->build().clone();
        if (bspline != nullptr)
        {
            functions.insert(bspline);
        }
    }
    return bspline;
}

void bspline_builder_delete(obj_ptr bspline_builder_ptr)
{
    auto func = get_bspline_builder(bspline_builder_ptr);

    if (func != nullptr)
    {
        bsplineBuilders.erase(bspline_builder_ptr);
        delete func;
    }
}

double *eval_row_major(obj_ptr function, double *x, int x_len)
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

double *eval_jacobian_row_major(obj_ptr function, double *x, int x_len)
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

double *eval_hessian_row_major(obj_ptr function, double *x, int x_len)
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

double *eval_col_major(obj_ptr function, double *x, int x_len)
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

            retVal = eval_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *eval_jacobian_col_major(obj_ptr function, double *x, int x_len)
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

            retVal = eval_jacobian_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }

    return retVal;
}

double *eval_hessian_col_major(obj_ptr function, double *x, int x_len)
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

            retVal = eval_hessian_row_major(func, row_major, x_len);
        }
        catch(const Exception &e)
        {
            set_error_string(e.what());
        }
        free(row_major);
    }
    return retVal;
}

int function_get_num_variables(obj_ptr function)
{
    int retVal = 0;

    auto func = get_function(function);
    if (func != nullptr)
    {
        retVal = func->getNumVariables();
    }

    return retVal;
}

void function_save(obj_ptr function, const char *filename)
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

void function_delete(obj_ptr function)
{
    auto func = get_function(function);

    if (func != nullptr)
    {
        functions.erase(function);
        delete func;
    }
}

}

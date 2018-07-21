/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_basis.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "kronecker_product.h"
#include <iostream>


namespace SPLINTER
{

BSplineBasis::BSplineBasis(std::vector<unsigned int> degrees, const std::vector< std::vector<double> > &knot_vectors)
    : num_variables(knot_vectors.size())
{
    if (knot_vectors.size() != degrees.size())
        throw Exception("BSplineBasis::BSplineBasis: Number of knot vectors is not equal to number of degrees.");

    // Set univariate bases
    bases.clear();
    for (unsigned int i = 0; i < num_variables; i++)
    {
        bases.emplace_back(BSplineBasis1D(degrees.at(i), knot_vectors.at(i)));

        // Adjust target number of basis functions used in e.g. refinement
        if (num_variables > 2)
        {
            // One extra knot is allowed
            bases.at(i).set_num_basis_functions_target((degrees.at(i) + 1) + 1); // Minimum degree+1
        }
    }
}

SparseVector BSplineBasis::eval(const DenseVector &x) const
{
    // Evaluate basis functions for each variable i and compute the tensor product of the function values
    std::vector<SparseVector> basis_function_values;

    for (int var = 0; var < x.size(); var++)
        basis_function_values.push_back(bases.at(var).eval(x(var)));

    return kronecker_product_vectors(basis_function_values);
}

// Old implementation of Jacobian
DenseMatrix BSplineBasis::eval_basis_jacobian_old(const DenseVector &x) const
{
    // Jacobian basis matrix
    DenseMatrix J;
    J.setZero(get_num_basis_functions(), num_variables);

    // Calculate partial derivatives
    for (unsigned int i = 0; i < num_variables; i++)
    {
        // One column in basis jacobian
        DenseVector bi; bi.setOnes(1);
        for (unsigned int j = 0; j < num_variables; j++)
        {
            DenseVector temp = bi;
            DenseVector xi;
            if (j == i)
            {
                // Differentiated basis
                xi = bases.at(j).eval_first_derivative(x(j));
            }
            else
            {
                // Normal basis
                xi = bases.at(j).eval(x(j));
            }

            bi = kroneckerProduct(temp, xi);
        }

        // Fill out column
        J.block(0,i,bi.rows(),1) = bi.block(0,0,bi.rows(),1);
    }

    return J;
}

// NOTE: does not pass tests
SparseMatrix BSplineBasis::eval_basis_jacobian(const DenseVector &x) const
{
    // Jacobian basis matrix
    SparseMatrix J(get_num_basis_functions(), num_variables);
    //J.setZero(num_basis_functions(), numInputs);

    // Calculate partial derivatives
    for (unsigned int i = 0; i < num_variables; ++i)
    {
        // One column in basis jacobian
        std::vector<SparseVector> values(num_variables);

        for (unsigned int j = 0; j < num_variables; ++j)
        {
            if (j == i)
            {
                // Differentiated basis
                values.at(j) = bases.at(j).eval_derivative(x(j), 1);
            }
            else
            {
                // Normal basis
                values.at(j) = bases.at(j).eval(x(j));
            }
        }

        SparseVector Ji = kronecker_product_vectors(values);

        // Fill out column i
        for (int k = 0; k < Ji.outerSize(); ++k)
        {
            for (SparseVector::InnerIterator it(Ji, k); it; ++it)
            {
                if (it.value() != 0)
                    J.insert(it.row(), i) = it.value();
            }
        }
        //J.block(0,i,Ji.rows(),1) = bi.block(0,0,Ji.rows(),1);
    }

    J.makeCompressed();

    return J;
}

SparseMatrix BSplineBasis::eval_basis_jacobian2(const DenseVector &x) const
{
    // Jacobian basis matrix
    SparseMatrix J(get_num_basis_functions(), num_variables);

    // Evaluate B-spline basis functions before looping
    std::vector<SparseVector> funcValues(num_variables);
    std::vector<SparseVector> gradValues(num_variables);

    for (unsigned int i = 0; i < num_variables; ++i)
    {
        funcValues[i] = bases.at(i).eval(x(i));
        gradValues[i] = bases.at(i).eval_first_derivative(x(i));
    }

    // Calculate partial derivatives
    for (unsigned int i = 0; i < num_variables; i++)
    {
        std::vector<SparseVector> values(num_variables);

        for (unsigned int j = 0; j < num_variables; j++)
        {
            if (j == i)
                values.at(j) = gradValues.at(j); // Differentiated basis
            else
                values.at(j) = funcValues.at(j); // Normal basis
        }

        SparseVector Ji = kronecker_product_vectors(values);

        // Fill out column
        for (SparseVector::InnerIterator it(Ji); it; ++it)
            J.insert(it.row(),i) = it.value();
    }

    return J;
}

SparseMatrix BSplineBasis::eval_basis_hessian(const DenseVector &x) const
{
    // Hessian basis matrix
    /* Hij = B1 x ... x DBi x ... x DBj x ... x Bn
     * (Hii = B1 x ... x DDBi x ... x Bn)
     * Where B are basis functions evaluated at x,
     * DB are the derivative of the basis functions,
     * and x is the kronecker product.
     * Hij is in R^(numBasisFunctions x 1)
     * so that basis hessian H is in R^(num_basis_functions*numInputs x numInputs)
     * The real B-spline Hessian is calculated as (c^T x 1^(numInputs x 1))*H
     */
    SparseMatrix H(get_num_basis_functions()*num_variables, num_variables);
    //H.setZero(num_basis_functions()*numInputs, numInputs);

    // Calculate partial derivatives
    // Utilizing that Hessian is symmetric
    // Filling out lower left triangular
    for (unsigned int i = 0; i < num_variables; i++) // row
    {
        for (unsigned int j = 0; j <= i; j++) // col
        {
            // One column in basis jacobian
            SparseMatrix Hi(1,1);
            Hi.insert(0,0) = 1;

            for (unsigned int k = 0; k < num_variables; k++)
            {
                SparseMatrix temp = Hi;
                SparseMatrix Bk;
                if (i == j && k == i)
                {
                    // Diagonal element
                    Bk = bases.at(k).eval_derivative(x(k), 2);
                }
                else if (k == i || k == j)
                {
                    Bk = bases.at(k).eval_derivative(x(k), 1);
                }
                else
                {
                    Bk = bases.at(k).eval(x(k));
                }
                Hi = kroneckerProduct(temp, Bk);
            }

            // Fill out column
            for (int k = 0; k < Hi.outerSize(); ++k)
            {
                for (SparseMatrix::InnerIterator it(Hi, k); it; ++it)
                {
                    if (it.value() != 0)
                    {
                        int row = i * get_num_basis_functions() + it.row();
                        int col = j;
                        H.insert(row, col) = it.value();
                    }
                }
            }
        }
    }

    H.makeCompressed();

    return H;
}

SparseMatrix BSplineBasis::insert_knots(double tau, unsigned int dim, unsigned int multiplicity)
{
    SparseMatrix A(1, 1);
    A.insert(0, 0) = 1;

    // Calculate multivariate knot insertion matrix
    for (unsigned int i = 0; i < num_variables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        if (i == dim)
        {
            // Build knot insertion matrix
            Ai = bases.at(i).insert_knots(tau, multiplicity);
        }
        else
        {
            // No insertion - identity matrix
            int m = bases.at(i).get_num_basis_functions();
            Ai.resize(m, m);
            Ai.setIdentity();
        }

//        A = kroneckerProduct(temp, Ai);
        A = my_kronecker_product(temp, Ai);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::refine_knots()
{
    SparseMatrix A(1, 1);
    A.insert(0, 0) = 1;

    for (unsigned int i = 0; i < num_variables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).refine_knots();

        //A = kroneckerProduct(temp, Ai);
        A = my_kronecker_product(temp, Ai);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::refine_knots_locally(const DenseVector &x)
{
    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < num_variables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).refine_knots_locally(x(i));

        //A = kroneckerProduct(temp, Ai);
        A = my_kronecker_product(temp, Ai);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::decompose_to_bezier()
{
    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < num_variables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).decompose_to_bezier();

        //A = kroneckerProduct(temp, Ai);
        A = my_kronecker_product(temp, Ai);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::reduce_support(const std::vector<double> &lb, const std::vector<double> &ub)
{
    if (lb.size() != ub.size() || lb.size() != num_variables)
        throw Exception("BSplineBasis::reduce_support: Incompatible dimension of domain bounds.");

    SparseMatrix A(1, 1);
    A.insert(0, 0) = 1;

    for (unsigned int i = 0; i < num_variables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).reduce_support(lb.at(i), ub.at(i));

        //A = kroneckerProduct(temp, Ai);
        A = my_kronecker_product(temp, Ai);
    }

    A.makeCompressed();

    return A;
}

std::vector<unsigned int> BSplineBasis::get_basis_degrees() const
{
    std::vector<unsigned int> degrees;
    for (const auto& basis : bases)
        degrees.push_back(basis.get_basis_degree());
    return degrees;
}

unsigned int BSplineBasis::get_basis_degree(unsigned int dim) const
{
    return bases.at(dim).get_basis_degree();
}

unsigned int BSplineBasis::get_num_basis_functions(unsigned int dim) const
{
    return bases.at(dim).get_num_basis_functions();
}

unsigned int BSplineBasis::get_num_basis_functions() const
{
    unsigned int prod = 1;
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        prod *= bases.at(dim).get_num_basis_functions();
    }
    return prod;
}

BSplineBasis1D BSplineBasis::get_single_basis(unsigned int dim)
{
    return bases.at(dim);
}

std::vector<double> BSplineBasis::get_knot_vector(int dim) const
{
    return bases.at(dim).get_knot_vector();
}

std::vector< std::vector<double> > BSplineBasis::get_knot_vectors() const
{
    std::vector< std::vector<double> > knots;
    for (unsigned int i = 0; i < num_variables; i++)
        knots.push_back(bases.at(i).get_knot_vector());
    return knots;
}

unsigned int BSplineBasis::get_knot_multiplicity(unsigned int dim, double tau) const
{
    return bases.at(dim).knot_multiplicity(tau);
}

std::vector<unsigned int> BSplineBasis::get_num_basis_functions_target() const
{
    std::vector<unsigned int> ret;
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        ret.push_back(bases.at(dim).get_num_basis_functions_target() );
    }
    return ret;
}

unsigned int BSplineBasis::num_supported() const
{
    unsigned int ret = 1;
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        ret *= (bases.at(dim).get_basis_degree() + 1);
    }
    return ret;
}

bool BSplineBasis::inside_support(const DenseVector &x) const
{
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        if (!bases.at(dim).is_supported(x(dim)))
        {
            return false;
        }
    }
    return true;
}

std::vector<double> BSplineBasis::get_support_lower_bound() const
{
    std::vector<double> lb;
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        std::vector<double> knots = bases.at(dim).get_knot_vector();
        lb.push_back(knots.front());
    }
    return lb;
}

std::vector<double> BSplineBasis::get_support_upper_bound() const
{
    std::vector<double> ub;
    for (unsigned int dim = 0; dim < num_variables; dim++)
    {
        std::vector<double> knots = bases.at(dim).get_knot_vector();
        ub.push_back(knots.back());
    }
    return ub;
}

} // namespace SPLINTER

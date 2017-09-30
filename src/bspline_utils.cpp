/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_utils.h"
#include <iostream>
#include <utilities.h>

namespace SPLINTER
{

SparseMatrix compute_basis_function_matrix(const BSpline &bspline, const DataTable &data)
{
    unsigned int num_samples = data.get_num_samples();

    // TODO: Reserve nnz per row (degree+1)
    //int nnzPrCol = bspline.basis.num_supported();

    SparseMatrix A(num_samples, bspline.get_num_basis_functions());
    //A.reserve(DenseVector::Constant(num_samples, nnzPrCol)); // TODO: should reserve nnz per row!

    int i = 0;
    for (auto it = data.cbegin(); it != data.cend(); ++it, ++i)
    {
        DenseVector xi = std_to_eig_vec(it->getX());
        SparseVector basis_values = bspline.eval_basis(xi);
        for (SparseVector::InnerIterator it2(basis_values); it2; ++it2)
            A.insert(i, it2.index()) = it2.value();
    }

    A.makeCompressed();

    return A;
}

DenseMatrix stack_sample_values(const DataTable &data)
{
    DenseMatrix B = DenseMatrix::Zero(data.get_num_samples(), data.get_dim_y());

    int i = 0;
    for (auto it = data.cbegin(); it != data.cend(); ++it, ++i)
    {
        auto y = it->getY();
        for (unsigned int j = 0; j < data.get_dim_y(); ++j)
            B(i, j) = y.at(j);
    }
    return B;
}

/**
 * Function for generating second order finite-difference matrix, which is used for penalizing the
 * (approximate) second derivative in control point calculation for P-splines.
 */
SparseMatrix compute_second_order_finite_difference_matrix(const BSpline &bspline)
{
    unsigned int num_variables = bspline.get_dim_x();

    // Number of (total) basis functions - defines the number of columns in D
    unsigned int num_cols = bspline.get_num_basis_functions();
    std::vector<unsigned int> num_basis_functions = bspline.get_num_basis_functions_per_variable();

    // Number of basis functions (and coefficients) in each variable
    std::vector<unsigned int> dims;
    for (unsigned int i = 0; i < num_variables; i++)
        dims.push_back(num_basis_functions.at(i));

    std::reverse(dims.begin(), dims.end());

    for (unsigned int i=0; i < num_variables; ++i)
        if (num_basis_functions.at(i) < 3)
            throw Exception("BSpline::Builder::getSecondOrderDifferenceMatrix: Need at least three coefficients/basis function per variable.");

    // Number of rows in D and in each block
    int num_rows = 0;
    std::vector< int > num_blk_rows;
    for (unsigned int i = 0; i < num_variables; i++)
    {
        int prod = 1;
        for (unsigned int j = 0; j < num_variables; j++)
        {
            if (i == j)
                prod *= (dims[j] - 2);
            else
                prod *= dims[j];
        }
        num_rows += prod;
        num_blk_rows.push_back(prod);
    }

    // Resize and initialize D
    SparseMatrix D(num_rows, num_cols);
    D.reserve(DenseVector::Constant(num_cols, 2*num_variables));   // D has no more than two elems per col per dim

    int i = 0; // Row index
    // Loop though each dimension (each dimension has its own block)
    for (unsigned int d = 0; d < num_variables; d++)
    {
        // Calculate left and right products
        int left_prod = 1;
        int right_prod = 1;
        for (unsigned int k = 0; k < d; k++)
        {
            left_prod *= dims[k];
        }
        for (unsigned int k = d+1; k < num_variables; k++)
        {
            right_prod *= dims[k];
        }

        // Loop through subblocks on the block diagonal
        for (int j = 0; j < right_prod; j++)
        {
            // Start column of current subblock
            int blkBaseCol = j*left_prod*dims[d];
            // Block rows [I -2I I] of subblock
            for (unsigned int l = 0; l < (dims[d] - 2); l++)
            {
                // Special case for first dimension
                if (d == 0)
                {
                    int k = j*left_prod*dims[d] + l;
                    D.insert(i,k) = 1;
                    k += left_prod;
                    D.insert(i,k) = -2;
                    k += left_prod;
                    D.insert(i,k) = 1;
                    i++;
                }
                else
                {
                    // Loop for identity matrix
                    for (int n = 0; n < left_prod; n++)
                    {
                        int k = blkBaseCol + l*left_prod + n;
                        D.insert(i,k) = 1;
                        k += left_prod;
                        D.insert(i,k) = -2;
                        k += left_prod;
                        D.insert(i,k) = 1;
                        i++;
                    }
                }
            }
        }
    }

    D.makeCompressed();

    return D;
}

SparseMatrix compute_weight_matrix(std::vector<double> weights)
{
    // TODO: use DiagonalMatrix here
    auto eig_weights = std_to_eig_vec(weights);
    DenseMatrix D = eig_weights.asDiagonal();
    return D.sparseView();
}

} // namespace SPLINTER
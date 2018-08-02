/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_utils.h"
#include <linear_solvers.h>
#include <iostream>
#include <utilities.h>

namespace SPLINTER
{

/**
 * Find coefficients of B-spline by solving:
 * min ||W*X*C - W*Y||^2 + alpha*||R||^2,
 * where
 * X = (m x n) matrix of n basis functions evaluated at m sample points,
 * Y = vector of m sample points y-values (or x-values when calculating knot averages),
 * C = B-spline coefficients (or knot averages),
 * R = Regularization matrix (n x n),
 * alpha = regularization parameter,
 * W = Diagonal weight matrix (m x m).
 *
 * The optimal control point matrix C is the solution of the linear system of equations:
 * (X'*W*X + alpha*R) C = X'*W*Y
 *
 * For performance gains, the various cases are solved as follows:
 * 1) No regularization or weighting: XC = Y
 * 2) Weighting, no regularization: WXC = WY
 * 3) Regularization, no weighting: (X'X + alpha*R)C = X'Y
 * 4) Regularization and weighting: (X'*W*X + alpha*R)C = X'*W*Y
 *
 * TODO: Use existing control points C0 as starting point and compute new control points as C = DeltaC + C0
 *
 */
DenseMatrix compute_control_points(const BSpline &bspline, const DataTable &data, BSpline::Smoothing smoothing,
                                   double alpha, std::vector<double> weights)
{
    unsigned int num_basis_functions = bspline.get_num_basis_functions();
    SparseMatrix X = compute_basis_function_matrix(bspline, data);
    DenseMatrix Y = stack_sample_values(data);

    SparseMatrix A; // Left-hand side matrix
    DenseMatrix B; // Right-hand side matrix

    // Weight matrix
    if (weights.size() > 0) {
        // NOTE: Consider using Eigen::DiagonalMatrix<double, num_samples> W()
        // SparseMatrix W(num_samples, num_samples);
        SparseMatrix W = compute_weight_matrix(weights);
        A = W*X;
        B = W*Y;
    } else {
        A = X;
        B = Y;
    }

    // Regularization
    if (smoothing != BSpline::Smoothing::NONE) {
        
        // Regularization matrix
        SparseMatrix R(num_basis_functions, num_basis_functions);

        if (smoothing == BSpline::Smoothing::IDENTITY) {
            /*
             * Tikhonov regularization (or ridge regression) with the Identity matrix
             * See: https://en.wikipedia.org/wiki/Tikhonov_regularization
             */
            auto I = SparseMatrix(num_basis_functions, num_basis_functions);
            I.setIdentity();
            R = I;
        }
        else if (smoothing == BSpline::Smoothing::PSPLINE) {
            /*
             * The P-Spline is a smoothing B-spline which relaxes the interpolation constraints on the control points to allow
             * smoother spline curves. It minimizes an objective which penalizes both deviation from sample points (to lower bias)
             * and the magnitude of second derivatives (to lower variance).
             *
             * Regularization matrix is given as R = D'*D, where D is the second-order finite difference matrix
             */
            SparseMatrix D = compute_second_order_finite_difference_matrix(bspline);
            R = D.transpose()*D;
        }
        
        // NOTE: Using eval to avoid aliasing
        // NOTE2: consider changing regularization factor to (alpha/numSample)
        A = (X.transpose()*A + alpha*R).eval();
        B = (X.transpose()*B).eval();
    }

    // Solve equation AC = B for control points C
    DenseMatrix C;

    // Solve as dense if both the number of rows and columns of A are less than 100
    bool solve_as_dense = (std::max(A.rows(), A.cols()) < 100);

    if (!solve_as_dense) {

        #ifndef NDEBUG
        std::cout << "compute_control_points: Computing B-spline control points using sparse solver." << std::endl;
        #endif // NDEBUG

        if (A.rows() == A.cols()) {
            // NOTE: If using SparseLU, the A matrix must be square
            SparseLU<DenseMatrix> s;
            solve_as_dense = !s.solve(A, B, C);
        }
        else {
            // Using SparseQR for non-square A matrix
            SparseQR<DenseMatrix> s;
            solve_as_dense = !s.solve(A, B, C);
        }
    }

    if (solve_as_dense) {

        #ifndef NDEBUG
        std::cout << "compute_control_points: Computing B-spline control points using dense solver." << std::endl;
        #endif // NDEBUG

        DenseMatrix Ad = A.toDense();
        DenseQR<DenseMatrix> s;
//        DenseSVD<DenseMatrix> s;
        if (!s.solve(Ad, B, C)) {
            throw Exception("compute_control_points: Failed to solve for B-spline control points.");
        }
    }

    return C;
}

SparseMatrix compute_basis_function_matrix(const BSpline &bspline, const DataTable &data)
{
    unsigned int num_samples = data.get_num_samples();

    SparseMatrix A(num_samples, bspline.get_num_basis_functions());

    // Reserve memory for A.
    // Assume ColMajor storage order, in which case the inner vectors of SparseMatrix are columns.
    // If the storage order ever changes, change the reserve statement to:
    //     A.reserve( Eigen::VectorXi::Constant( num_samples, bspline.get_num_supported() ) );
    static_assert( A.IsRowMajor == false, "" );
    // Although the number of non-zero elements per row equals bspline.get_num_supported(),
    // SparseMatrix::reserve requires a vector of size SparseMatrix::cols() a an argument;
    Eigen::VectorXi reserve_sizes( A.cols() );
    // Assuming that non-zero elements in A are distributed randomly, the average number of non-zero elements per column is:
    const double average_nnz_per_col = num_samples * bspline.get_num_supported() / bspline.get_num_basis_functions();
    // Add safety margin of sqrt(N), but don't exceed num_samples ( which otherwise would happen if e.g. bspline.get_num_supported() == bspline.get_num_basis_functions() )
    const int reserve_col_size = std::min( static_cast<int>( average_nnz_per_col + sqrt( average_nnz_per_col ) ), static_cast<int>(num_samples) );
    reserve_sizes.fill( reserve_col_size );
    A.reserve( reserve_sizes );

    int i = 0;
    for (auto it = data.cbegin(); it != data.cend(); ++it, ++i)
    {
        DenseVector xi = std_to_eig_vec(it->get_x());
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
        auto y = it->get_y();
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
            throw Exception("compute_second_order_finite_difference_matrix: Need at least three coefficients/basis functions per variable.");

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

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_builder.h"
#include "bspline_utils.h"
#include <linear_solvers.h>
#include <serializer.h>
#include <iostream>
#include <utilities.h>

namespace SPLINTER
{
// Default constructor
BSpline::Builder::Builder(unsigned int dim_x, unsigned int dim_y)
        : _dim_x(dim_x),
          _dim_y(dim_y),
          _degrees(std::vector<unsigned int>(_dim_x, 3)),
          _numBasisFunctions(std::vector<unsigned int>(_dim_x, 1)),
          _knotSpacing(KnotSpacing::AS_SAMPLED)
{
}

/**
 * Fit B-spline to data
 */
BSpline BSpline::Builder::fit(const DataTable &data,
                              Smoothing smoothing,
                              double alpha,
                              std::vector<double> weights) const
{
    if (data.getDimX() != _dim_x)
        throw Exception("BSpline::Builder::fit: Expected " + std::to_string(_dim_x) + " input variables.");

    if (data.getDimY() != _dim_y)
        throw Exception("BSpline::Builder::fit: Expected " + std::to_string(_dim_y) + " output variables.");

    if (alpha < 0)
        throw Exception("BSpline::Builder::fit: alpha must be non-negative.");

    if (weights.size() > 0 && data.getNumSamples() != weights.size()) {
        throw Exception("BSpline::Builder::fit: number of weights must equal number of data points.");
    }

#ifndef NDEBUG
    if (!data.isGridComplete())
        std::cout << "BSpline::Builder::fit: Building B-spline from irregular (incomplete) grid." << std::endl;
#endif // NDEBUG

    // Build knot vectors
    auto knotVectors = computeKnotVectors(data, _degrees, _numBasisFunctions, _knotSpacing);

    // Build B-spline (with default coefficients)
    auto bspline = BSpline(_dim_x, _dim_y, knotVectors, _degrees);

    // Compute coefficients from samples and update B-spline
    auto coefficients = computeControlPoints(bspline, data, smoothing, alpha, weights);
    bspline.setControlPoints(coefficients);

    return bspline;
}

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
 */
DenseMatrix BSpline::Builder::computeControlPoints(const BSpline &bspline,
                                                   const DataTable &data,
                                                   Smoothing smoothing,
                                                   double alpha,
                                                   std::vector<double> weights) const
{
    unsigned int num_samples = data.getNumSamples();
    unsigned int num_basis_functions = bspline.getNumBasisFunctions();
    SparseMatrix X = computeBasisFunctionMatrix(bspline, data);
    SparseMatrix Xt = X.transpose();
    DenseMatrix Y = stackSamplePointValues(data);

    // Regularization matrix
    SparseMatrix R(num_basis_functions, num_basis_functions);

    if (smoothing == Smoothing::IDENTITY) {
        /*
         * Tikhonov regularization (or ridge regression) with the Identity matrix
         * See: https://en.wikipedia.org/wiki/Tikhonov_regularization
         */
        auto I = SparseMatrix(num_basis_functions, num_basis_functions);
        I.setIdentity();
        R = I;
    }
    else if (smoothing == Smoothing::PSPLINE)
    {
        /*
         * The P-Spline is a smoothing B-spline which relaxes the interpolation constraints on the control points to allow
         * smoother spline curves. It minimizes an objective which penalizes both deviation from sample points (to lower bias)
         * and the magnitude of second derivatives (to lower variance).
         *
         * Regularization matrix is given as R = D'*D, where D is the second-order finite difference matrix
         */

        // Second order finite difference matrix
        SparseMatrix D = computeSecondOrderFiniteDifferenceMatrix(bspline);
        R = D.transpose()*D;
    }

    // Weight matrix
    // NOTE: Consider using Eigen::DiagonalMatrix<double, num_samples> W()
    SparseMatrix W(num_samples, num_samples);

    if (weights.size() > 0) {
        W = computeWeightMatrix(weights);
    } else {
        W.setIdentity();
    }

    // Left-hand side matrix
    // NOTE2: consider changing regularization factor to (alpha/numSample)
    SparseMatrix A = Xt*W*X;

    if (smoothing != Smoothing::NONE) {
        A += alpha*R;
    }

    // Compute right-hand side matrices
    DenseMatrix B = Xt*W*Y;

    // Solve equation AC = B for control points C
    DenseMatrix C;

    int num_equations = A.rows();
    int max_num_equations = 100;
    bool solve_as_dense = (num_equations < max_num_equations);

    if (!solve_as_dense)
    {
        #ifndef NDEBUG
        std::cout << "BSpline::Builder::computeBSplineCoefficients: Computing B-spline control points using sparse solver." << std::endl;
        #endif // NDEBUG

        SparseLU<DenseMatrix> s;

        solve_as_dense = !s.solve(A, B, C);
    }

    if (solve_as_dense)
    {
        #ifndef NDEBUG
        std::cout << "BSpline::Builder::computeBSplineCoefficients: Computing B-spline control points using dense solver." << std::endl;
        #endif // NDEBUG

        DenseMatrix Ad = A.toDense();
        DenseQR<DenseMatrix> s;
        // DenseSVD<DenseVector> s;
        if (!s.solve(Ad, B, C))
        {
            throw Exception("BSpline::Builder::computeBSplineCoefficients: Failed to solve for B-spline coefficients.");
        }
    }

    return C;
}

} // namespace SPLINTER
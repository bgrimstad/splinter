/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_UTILS_H
#define SPLINTER_BSPLINE_UTILS_H

#include <data_table.h>
#include <bspline.h>

namespace SPLINTER
{

// Control point computations
DenseMatrix compute_control_points(const BSpline &bspline, const DataTable &data, BSpline::Smoothing smoothing,
                                   double alpha, std::vector<double> weights);

// Matrix of basis functions evaluated at samples
SparseMatrix compute_basis_function_matrix(const BSpline &bspline, const DataTable &data);

// Stack samples in DenseMatrix
DenseMatrix stack_sample_values(const DataTable &data);

// P-spline control point calculation
SparseMatrix compute_second_order_finite_difference_matrix(const BSpline &bspline);

// Compute weights matrix from weight vector
SparseMatrix compute_weight_matrix(std::vector<double> weights);

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_UTILS_H

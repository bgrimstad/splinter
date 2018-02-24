/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_BUILDER_H
#define SPLINTER_BSPLINE_BUILDER_H

#include <data_table.h>
#include <knot_builders.h>
#include <bspline.h>

namespace SPLINTER
{

/**
 * Convenience functions for B-spline fitting
 */


/**
 * Create a B-spline that interpolates the sample points.
 * @param data A table of sample points on a regular grid.
 * @param degree The degree of the B-spline basis functions. Default degree is 3 (cubic).
 * @return A B-spline that interpolates the sample points.
 */
BSpline bspline_interpolator(const DataTable &data, unsigned int degree = 3);


/**
 * Create a B-spline that smooths the sample points using regularization (weight decay).
 * @param data A table of sample points on a regular grid.
 * @param degree The degree of the B-spline basis functions. Default degree is 3 (cubic).
 * @param smoothing Type of regularization to use - see BSpline::Smoothing. Default is smoothing is BSpline::PSLINE.
 * @param alpha Smoothing/regularization factor.
 * @param weights Sample weights.
 * @return A B-spline that smooths the sample points.
 */
BSpline bspline_smoother(const DataTable &data, unsigned int degree = 3,
                         BSpline::Smoothing smoothing = BSpline::Smoothing::PSPLINE, double alpha = 0.1,
                         std::vector<double> weights = std::vector<double>());


/**
 * Create a unfitted (zero-valued) B-spline. This builder gives the user more control and allows for construction of
 * B-splines with non-default knot vectors and different degrees for the univariate basis functions.
 * @param data A table of sample points on a regular grid.
 * @param degree The degrees of the B-spline basis functions.
 * @param knot_spacing The knot spacing method used to build knot vectors.
 * @param num_basis_functions The desired number of basis functions in each univariate basis (must be at least degree + 1).
 * @return A B-spline that smooths the sample points.
 */
BSpline bspline_unfitted(const DataTable &data, const std::vector<unsigned int> &degrees, KnotSpacing knot_spacing,
                         const std::vector<unsigned int> &num_basis_functions);


} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BUILDER_H

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
#include <serializer.h>

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
BSpline BSpline::Builder::fit(const DataTable &data, Smoothing smoothing, double alpha,
                              std::vector<double> weights) const
{
    if (data.get_dim_x() != _dim_x)
        throw Exception("BSpline::Builder::fit: Expected " + std::to_string(_dim_x) + " input variables.");

    if (data.get_dim_y() != _dim_y)
        throw Exception("BSpline::Builder::fit: Expected " + std::to_string(_dim_y) + " output variables.");

    if (alpha < 0)
        throw Exception("BSpline::Builder::fit: alpha must be non-negative.");

    if (weights.size() > 0 && data.get_num_samples() != weights.size()) {
        throw Exception("BSpline::Builder::fit: number of weights must equal number of data points.");
    }

    // Build knot vectors
    auto knotVectors = compute_knot_vectors(data, _degrees, _knotSpacing, _numBasisFunctions);

    // Build B-spline (with default coefficients)
    return BSpline(_dim_x, _dim_y, knotVectors, _degrees).fit(data, smoothing, alpha, weights);
}

/**
 * Create a B-spline that interpolates the sample points
 * @param data A table of sample points on a regular grid
 * @param degree The degree of the B-spline basis functions
 * @return A B-spline that interpolates the sample points
 */
BSpline bspline_interpolator(const DataTable &data, unsigned int degree)
{
    auto dim_x = data.get_dim_x();
    auto dim_y = data.get_dim_y();
    auto degrees = std::vector<unsigned int>(dim_x, degree);
    auto knot_spacing = KnotSpacing::AS_SAMPLED;
    auto knot_vectors = compute_knot_vectors(data, degrees, knot_spacing);

    return BSpline(dim_x, dim_y, knot_vectors, degrees).fit(data);
}

/**
 * Create a cubic B-spline that interpolates the sample points
 * In the multivariate case, the multivariate basis functions will be products of univariate cubic basis functions,
 * i.e. for two variables the basis functions are bi-cubic, etc.
 * @param data A table of sample points on a regular grid
 * @return A cubic B-spline that interpolates the sample points
 */
BSpline cubic_bspline_interpolator(const DataTable &data)
{
    return bspline_interpolator(data, 3);
}

/**
 * Create a P-spline (penalized B-spline) that smooths the sample points
 * @param data A table of sample points on a regular grid
 * @param degree The degree of the B-spline basis functions
 * @param alpha Smoothing/regularization factor
 * @return A B-spline that smooths the sample points (P-spline)
 */
BSpline pspline_approximator(const DataTable &data, unsigned int degree, double alpha)
{
    auto dim_x = data.get_dim_x();
    auto dim_y = data.get_dim_y();
    auto degrees = std::vector<unsigned int>(dim_x, degree);
    auto knot_spacing = KnotSpacing::AS_SAMPLED;
    auto knot_vectors = compute_knot_vectors(data, degrees, knot_spacing);

    return BSpline(dim_x, dim_y, knot_vectors, degrees).fit(data, BSpline::Smoothing::PSPLINE, alpha);
}

/**
 * Create a cubic P-spline (penalized B-spline) that smooths the sample points.
 * In the multivariate case, the multivariate basis functions will be products of univariate cubic basis functions,
 * i.e. for two variables the basis functions are bi-cubic, etc.
 * @param data A table of sample points on a regular grid
 * @param alpha Smoothing/regularization factor
 * @return A cubic B-spline that smooths the sample points (cubic P-spline)
 */
BSpline cubic_pspline_approximator(const DataTable &data, double alpha)
{
    return pspline_approximator(data, 3, alpha);
}


} // namespace SPLINTER
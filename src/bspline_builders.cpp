/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_builders.h"


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
    auto knotVectors = build_knot_vectors(data, _degrees, _knotSpacing, _numBasisFunctions);

    // Build B-spline (with default coefficients)
    return BSpline(_degrees, knotVectors, _dim_y).fit(data, smoothing, alpha, weights);
}

BSpline bspline_interpolator(const DataTable &data, unsigned int degree)
{
    auto dim_x = data.get_dim_x();
    auto dim_y = data.get_dim_y();
    auto degrees = std::vector<unsigned int>(dim_x, degree);
    auto knot_vectors = build_knot_vectors(data, degrees);

    return BSpline(degrees, knot_vectors, dim_y).fit(data);
}

BSpline bspline_smoother(const DataTable &data, unsigned int degree, BSpline::Smoothing smoothing, double alpha,
                         std::vector<double> weights)
{
    auto dim_x = data.get_dim_x();
    auto dim_y = data.get_dim_y();
    auto degrees = std::vector<unsigned int>(dim_x, degree);
    auto knot_vectors = build_knot_vectors(data, degrees);

    return BSpline(degrees, knot_vectors, dim_y).fit(data, smoothing, alpha, weights);
}

BSpline bspline_unfitted(const DataTable &data, const std::vector<unsigned int> &degrees, KnotSpacing knot_spacing,
                         const std::vector<unsigned int> &num_basis_functions)
{
    auto dim_x = data.get_dim_x();
    auto dim_y = data.get_dim_y();
    if (dim_x != num_basis_functions.size())
        throw Exception("Size of num_basis_functions " + std::to_string(num_basis_functions.size())
                        + " does not match the input dimension " + std::to_string(dim_x) + " of the samples.");

    if (dim_y != degrees.size())
        throw Exception("Size of degrees " + std::to_string(degrees.size()) + " does not match the output dimensions "
                        + std::to_string(dim_y) + " of the samples.");

    auto knot_vectors = build_knot_vectors(data, degrees, knot_spacing, num_basis_functions);
    return BSpline(degrees, knot_vectors, dim_y);
}

} // namespace SPLINTER
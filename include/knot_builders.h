/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_KNOT_UTILS_H
#define SPLINTER_KNOT_UTILS_H

#include <vector>
#include <data_table.h>

namespace SPLINTER
{

/**
 * B-spline knot spacing
 */
enum class KnotSpacing {
    /*
     * Experimental knot spacing for testing purposes only.
     * Currently, it computes equidistant knots on an expanded interval.
     * NOTE: may change in the future
     */
    EXPERIMENTAL,

    /*
     * Clamped and with knots that mimic the spacing of sample points using a moving average filter.
     * Resulting knot vectors have p+1 multiplicity of end knots. This is the default knot spacing for B-spline
     * interpolation and smoothing.
     * Example for samples at 0, 1, ..., 5:
     *      Then, for degree 3 the knot vector becomes: [0, 0, 0, 0, 2, 3, 5, 5, 5, 5]
     *      and for degree 1 the knot vector becomes: [0, 0, 1, 2, 3, 4, 5, 5]
     */
    AS_SAMPLED,

    /*
     * Clamped knot vector with equidistant internal knots. p+1 multiplicity of end knots.
     * Example for samples at 0, 1, ..., 5:
     *      For degree 3: [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
     *      For degree 1: [0, 0, 1, 2, 3, 4, 5, 5]
     */
    EQUIDISTANT_CLAMPED,

    /*
     * Simple knot vector with equidistant knots. Single multiplicity of end knots. Does not vary with degree.
     * Example for samples at 0, 1, ..., 5:
     *      For degree 3: [0, 1, 2, 3, 4, 5]
     *      For degree 1: [0, 1, 2, 3, 4, 5]
     */
    EQUIDISTANT
};

/**
 * Computing knot vectors
 */
std::vector<std::vector<double>> build_knot_vectors(const DataTable &data, const std::vector<unsigned int> &degrees);

std::vector<std::vector<double>> build_knot_vectors(const DataTable &data, const std::vector<unsigned int> &degrees,
                                                    KnotSpacing knot_spacing,
                                                    const std::vector<unsigned int> &num_basis_functions);

/**
 * Computing single knot vector
 */
std::vector<double> build_knot_vector(const std::vector<double> &values, unsigned int degree,
                                      unsigned int num_basis_functions, KnotSpacing knot_spacing);

/**
 * Construct knot vector using moving average filter (AS_SAMPLED)
 */
std::vector<double> knot_vector_moving_average(const std::vector<double> &values, unsigned int degree);

/**
 * Compute clamped, equidistant knot vector (EQUIDISTANT_CLAMPED)
 */
std::vector<double> knot_vector_equidistant_clamped(const std::vector<double> &values, unsigned int degree,
                                                    unsigned int num_basis_functions = 0);

/**
 * Construct equidistant knot vector (EQUIDISTANT)
 */
std::vector<double> knot_vector_equidistant(const std::vector<double> &values, unsigned int degree,
                                            unsigned int num_basis_functions);

/**
 * Construct expanded equidistant knot vector (EXPERIMENTAL)
 */
std::vector<double> knot_vector_expanded_equidistant(const std::vector<double> &values, unsigned int degree,
                                                     unsigned int num_basis_functions);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_UTILS_H

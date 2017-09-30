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
     * Clamped and with knots that mimic the spacing of sample points using a moving average filter.
     * p+1 multiplicity of end knots.
     * Example:
     *      Given sample points [0, 1, 2, 3, 4, 5]
     *      Then, for degree 3 the knot vector becomes: [0, 0, 0, 0, 2, 3, 5, 5, 5, 5]
     *      and for degree 1 the knot vector becomes: [0, 0, 1, 2, 3, 4, 5, 5]
     */
            AS_SAMPLED,

    /*
     * Clamped knot vector with equidistant internal knots. p+1 multiplicity of end knots.
     * Example:
     *      Given samples on the interval [0, 5]
     *      For degree 3: [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
     *      For degree 1: [0, 0, 1, 2, 3, 4, 5, 5]
     */
            EQUIDISTANT,

    /*
     * Experimental knot spacing for testing purposes only.
     * Currently, it gives a non-clamped knot vector of equidistant knots.
     * NOTE: may eventually become EQUIDISTANT_NOT_CLAMPED
     * Example:
     *      For any degree: [0, 1, 2, 3, 4, 5]
     *
     */
            EXPERIMENTAL
};

/**
 * Computing knot vectors
 */
std::vector<std::vector<double>> compute_knot_vectors(const DataTable &data, std::vector<unsigned int> degrees,
                                                      KnotSpacing knot_spacing);

std::vector<std::vector<double>> compute_knot_vectors(const DataTable &data, std::vector<unsigned int> degrees,
                                                      KnotSpacing knot_spacing,
                                                      std::vector<unsigned int> num_basis_functions);
/**
 * Computing single knot vector
 */
std::vector<double> compute_knot_vector(const std::vector<double> &values, unsigned int degree,
                                        unsigned int num_basis_functions, KnotSpacing knot_spacing);

/**
 * Construct knot vector using moving average filter (AS_SAMPLED)
 */
std::vector<double> knot_vector_moving_average(const std::vector<double> &values, unsigned int degree);

/**
 * Compute clamped, equidistant knot vector (EQUIDISTANT)
 */
std::vector<double> knot_vector_equidistant(const std::vector<double> &values, unsigned int degree,
                                            unsigned int num_basis_functions = 0);

/**
 * Construct equidistant knot vector (EXPERIMENTAL)
 */
std::vector<double> knot_vector_equidistant_not_clamped(const std::vector<double> &values, unsigned int degree,
                                                        unsigned int num_basis_functions);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_UTILS_H

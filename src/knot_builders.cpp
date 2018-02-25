/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <knot_builders.h>
#include <utilities.h>
#include <algorithm>


namespace SPLINTER
{

// Compute all knot vectors from sample data
std::vector<std::vector<double>> build_knot_vectors(const DataTable &data, const std::vector<unsigned int> &degrees)
{
    auto num_basis_functions = std::vector<unsigned int>(data.get_dim_x(), 10); // Not used when knot spacing is AS_SAMPLED
    return build_knot_vectors(data, degrees, KnotSpacing::AS_SAMPLED, num_basis_functions);
}

// Compute all knot vectors from sample data
std::vector<std::vector<double>> build_knot_vectors(const DataTable &data,
                                                    const std::vector<unsigned int> &degrees,
                                                    KnotSpacing knot_spacing,
                                                    const std::vector<unsigned int> &num_basis_functions)
{
    auto dim_x = data.get_dim_x();

    if (dim_x != degrees.size() || dim_x != num_basis_functions.size())
        throw Exception("BSpline::Builder::build_knot_vectors: Inconsistent sizes on input vectors.");

    std::vector<std::vector<double>> grid = data.get_table_x();

    std::vector<std::vector<double>> knot_vectors;

    for (unsigned int i = 0; i < dim_x; ++i)
    {
        // Compute knot vector
        auto knot_vector = build_knot_vector(grid.at(i), degrees.at(i), num_basis_functions.at(i), knot_spacing);
        knot_vectors.push_back(knot_vector);
    }

    return knot_vectors;
}

// Compute a single knot vector from sample grid and degree
std::vector<double> build_knot_vector(const std::vector<double> &values, unsigned int degree,
                                      unsigned int num_basis_functions, KnotSpacing knot_spacing)
{
    switch (knot_spacing)
    {
        case KnotSpacing::EXPERIMENTAL:
            return knot_vector_expanded_equidistant(values, degree, num_basis_functions);
        case KnotSpacing::AS_SAMPLED:
            return knot_vector_moving_average(values, degree);
        case KnotSpacing::EQUIDISTANT_CLAMPED:
            return knot_vector_equidistant_clamped(values, degree, num_basis_functions);
        case KnotSpacing::EQUIDISTANT:
            return knot_vector_equidistant(values, degree, num_basis_functions);
        default:
            return knot_vector_moving_average(values, degree);
    }
}

/*
* Automatic construction of (p+1)-regular knot vector using moving average.
*
* Requirement:
* Knot vector should be of size n+p+1.
* End knots are should be repeated p+1 times.
*
* Computed sizes:
* n+2*(p) = n + p + 1 + (p - 1)
* k = (p - 1) values must be removed from sample vector.
* w = k + 3 window size in moving average
*
* Algorithm:
* 1) compute n - k values using moving average with window size w
* 2) repeat first and last value p + 1 times
*
* The resulting knot vector has n - k + 2*p = n + p + 1 knots.
*
* NOTE:
* For equidistant samples, the resulting knot vector is identically to
* the free end conditions knot vector used in cubic interpolation.
* That is, samples (a,b,c,d,e,f) produces the knot vector (a,a,a,a,c,d,f,f,f,f) for p = 3.
* For p = 1, (a,b,c,d,e,f) becomes (a,a,b,c,d,e,f,f).
*
* TODO:
* Does not work well when number of knots is << number of samples! For such cases
* almost all knots will lie close to the left samples. Try a bucket approach, where the
* samples are added to buckets and the knots computed as the average of these.
*/
std::vector<double> knot_vector_moving_average(const std::vector<double> &values, unsigned int degree)
{
    // Sort and remove duplicates
    std::vector<double> unique = extract_unique_sorted(values);

    // Compute sizes
    auto n = static_cast<unsigned int>(unique.size());
    auto k = degree-1; // knots to remove
    auto w = k + 3; // Window size

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knot_vector_moving_average: Only " << n
        << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
        << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    std::vector<double> knots(n-k-2, 0);

    // Compute (n-k-2) interior knots using moving average
    for (unsigned int i = 0; i < n-k-2; ++i)
    {
        double ma = 0;
        for (unsigned int j = 0; j < w; ++j)
            ma += unique.at(i+j);

        knots.at(i) = ma/w;
    }

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (auto i = 0u; i < degree + 1; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (auto i = 0u; i < degree + 1; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == unique.size() + degree + 1);

    return knots;
}

std::vector<double> knot_vector_equidistant_clamped(const std::vector<double> &values, unsigned int degree,
                                                    unsigned int num_basis_functions)
{
    // Sort and remove duplicates
    std::vector<double> unique = extract_unique_sorted(values);

    // Compute sizes
    auto n = static_cast<unsigned int>(unique.size());
    if (num_basis_functions > 0)
        n = num_basis_functions;
    auto k = degree-1; // knots to remove

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knot_vector_equidistant_clamped: Only " << n
        << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
        << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    // Compute (n-k-2) equidistant interior knots
    auto num_internal_knots = std::max(n-k-2, (unsigned int)0);
    num_internal_knots = std::min((unsigned int)10, num_internal_knots);
    std::vector<double> knots = linspace(unique.front(), unique.back(), num_internal_knots);

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (auto i = 0u; i < degree; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (auto i = 0u; i < degree; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

std::vector<double> knot_vector_equidistant(const std::vector<double> &values, unsigned int degree,
                                            unsigned int num_basis_functions)
{
    // Sort and remove duplicates
    std::vector<double> unique = extract_unique_sorted(values);

    // Build equidistant knots
    unsigned int nk = num_basis_functions + degree + 1;
    return linspace(unique.front(), unique.back(), nk);
}

std::vector<double> knot_vector_expanded_equidistant(const std::vector<double> &values, unsigned int degree,
                                                     unsigned int num_basis_functions)
{
    // Sort and remove duplicates
    std::vector<double> unique = extract_unique_sorted(values);

    // Number of knots
    unsigned int nk = num_basis_functions + degree + 1;

    // Expand each boundary by 10% of sample interval
    double lb = unique.front();
    double ub = unique.back();
    double delta = ub - lb;
    double expansion = 0.1; // A non-negative number
    return linspace(unique.front() - expansion*delta, unique.back() + expansion*delta, nk);
}

} // namespace SPLINTER

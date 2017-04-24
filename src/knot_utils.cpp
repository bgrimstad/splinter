/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <knot_utils.h>
#include <utilities.h>
#include <algorithm>

namespace SPLINTER
{

// Compute all knot vectors from sample data
std::vector<std::vector<double>> computeKnotVectors(const DataTable &data,
                                                    std::vector<unsigned int> degrees,
                                                    std::vector<unsigned int> num_basis_functions,
                                                    KnotSpacing knot_spacing)
{
    auto dim_x = data.getDimX();

    if (dim_x != degrees.size() || dim_x != num_basis_functions.size())
        throw Exception("BSpline::Builder::computeKnotVectors: Inconsistent sizes on input vectors.");

    std::vector<std::vector<double>> grid = data.getTableX();

    std::vector<std::vector<double>> knotVectors;

    for (unsigned int i = 0; i < dim_x; ++i)
    {
        // Compute knot vector
        auto knotVec = computeKnotVector(grid.at(i), degrees.at(i), num_basis_functions.at(i), knot_spacing);

        knotVectors.push_back(knotVec);
    }

    return knotVectors;
}

// Compute a single knot vector from sample grid and degree
std::vector<double> computeKnotVector(const std::vector<double> &values,
                                      unsigned int degree,
                                      unsigned int num_basis_functions,
                                      KnotSpacing knot_spacing)
{
    switch (knot_spacing)
    {
        case KnotSpacing::AS_SAMPLED:
            return knotVectorMovingAverage(values, degree);
        case KnotSpacing::EQUIDISTANT:
            return knotVectorEquidistant(values, degree, num_basis_functions);
        case KnotSpacing::EXPERIMENTAL:
            return knotVectorEquidistantNotClamped(values, degree, num_basis_functions);
        default:
            return knotVectorMovingAverage(values, degree);
    }
}

std::vector<double> knotVectorEquidistantNotClamped(const std::vector<double> &values,
                                                    unsigned int degree,
                                                    unsigned int numBasisFunctions)
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // Number of knots
    unsigned int nk = numBasisFunctions + degree + 1;

    // Compute (n-k-2) equidistant interior knots
    double lb = unique.front();
    double ub = unique.back();
    double delta = ub - lb;
    double expansion = 0.1; // A non-negative number
    std::vector<double> knots = linspace(unique.front() - expansion*delta, unique.back() + expansion*delta, nk);

    return knots;
}

/*
* Automatic construction of (p+1)-regular knot vector
* using moving average.
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
std::vector<double> knotVectorMovingAverage(const std::vector<double> &values,
                                            unsigned int degree)
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // Compute sizes
    unsigned int n = unique.size();
    unsigned int k = degree-1; // knots to remove
    unsigned int w = k + 3; // Window size

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knotVectorMovingAverage: Only " << n
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
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

std::vector<double> knotVectorEquidistant(const std::vector<double> &values,
                                          unsigned int degree,
                                          unsigned int numBasisFunctions = 0)
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // Compute sizes
    unsigned int n = unique.size();
    if (numBasisFunctions > 0)
        n = numBasisFunctions;
    unsigned int k = degree-1; // knots to remove

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knotVectorMovingAverage: Only " << n
        << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
        << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    // Compute (n-k-2) equidistant interior knots
    unsigned int numIntKnots = std::max(n-k-2, (unsigned int)0);
    numIntKnots = std::min((unsigned int)10, numIntKnots);
    std::vector<double> knots = linspace(unique.front(), unique.back(), numIntKnots);

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (unsigned int i = 0; i < degree; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

} // namespace SPLINTER

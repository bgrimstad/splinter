/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "test_utils.h"
#include <utilities.h>
#include <utils/test_function_collection.h>
#include <random>
#include <iostream>
#include <Catch.h>


using namespace std;

namespace SPLINTER
{

// See https://en.wikipedia.org/wiki/Relative_change_and_difference#Formulae
double rel_error(double exactVal, double approxVal)
{
    double maxAbsVal = std::max(std::abs(exactVal), std::abs(approxVal));

    // Both are ~0
    if (maxAbsVal < 1e-14)
    {
        return 0.0;
    }

    double absError = std::abs(exactVal - approxVal);

    return std::min(absError, absError / maxAbsVal);
}

std::vector<double> get_random_vector(unsigned int length, unsigned int seed)
{
    // Random number generator (using normal distribution)
    std::default_random_engine generator(seed); // Seeded to be reproducible
    double mean = 0.;
    double stddev = 10.;
    std::normal_distribution<double> distribution(mean, stddev);

    std::vector<double> random_vector;
    for (unsigned int i = 0; i < length; ++i) {
        random_vector.push_back(distribution(generator));
    }

    return random_vector;
}

double six_hump_camel_back(std::vector<double> x)
{
    assert(x.size() == 2);
    double x0 = x.at(0);
    double x1 = x.at(1);
    return (4 - 2.1*x0*x0 + (1/3.)*x0*x0*x0*x0)*x0*x0 + x0*x1 + (-4 + 4*x1*x1)*x1*x1;
}

double matrix_one_norm(const DenseMatrix &mat)
{
    return mat.lpNorm<1>();
}

double matrix_two_norm(const DenseMatrix &mat)
{
    return mat.lpNorm<2>();
}

double matrix_inf_norm(const DenseMatrix &mat)
{
    return mat.lpNorm<Eigen::Infinity>();
}

double log(double base, double x)
{
    return std::log(x) / std::log(base);
}

// Enumerates all permutations
std::vector<std::vector<double>> multi_linspace(std::vector<double> start, std::vector<double> end,
                                                std::vector<unsigned int> points)
{
    assert(start.size() == end.size() && end.size() == points.size());

    size_t nDims = start.size();
    size_t numPoints = 1;
    for(size_t i = 0; i < nDims; i++) {
        numPoints *= points.at(i);
    }

#ifndef NDEBUG
    if(numPoints > 10000) {
        std::cout << "Warning: Enumerating " << numPoints << " points." << std::endl;
    }
#endif // ifndef NDEBUG

    auto result = std::vector<std::vector<double>>(numPoints);

    auto dx = std::vector<double>(nDims);
    auto it = std::vector<unsigned int>(nDims);
    auto cur = std::vector<double>(nDims);

    for(unsigned int i = 0; i < nDims; i++) {
        cur.at(i) = start.at(i);
        it.at(i) = 0;
        dx.at(i) = (end.at(i) - start.at(i)) / (points.at(i) - 1);
    }

    size_t curDim = 0;
    size_t i = 0;
    // Add the start vector
    result.at(i++) = std::vector<double>(start);
    while(true) {
        curDim = 0;
        while(curDim < nDims && it.at(curDim)+1 >= points.at(curDim)) {
            it.at(curDim) = 0;
            cur.at(curDim) = start.at(curDim) + dx.at(curDim) * it.at(curDim);
            curDim++;
        }
        // If we're unable to find a number that can be increased,
        // it means were at the highest number
        // (for 4 digits in decimal that is 9999)
        if(curDim >= nDims) {
            break;
        }

        it.at(curDim)++;
        cur.at(curDim) = start.at(curDim) + dx.at(curDim) * it.at(curDim);

        result.at(i++) = std::vector<double>(cur);
    }

    assert(i == numPoints);

    return result;
}

std::vector<std::vector<double>> multi_linspace(int dim, double start, double end, unsigned int points)
{
    auto startVec = std::vector<double>(dim);
    auto endVec = std::vector<double>(dim);
    auto pointsVec = std::vector<unsigned int>(dim);

    for(int i = 0; i < dim; i++) {
        startVec.at(i) = start;
        endVec.at(i) = end;
        pointsVec.at(i) = points;
    }

    return multi_linspace(startVec, endVec, pointsVec);
}

std::vector<std::vector<double>> multi_linspace(int dim)
{
    // Total number of points to evaluate
    // Will be distributed evenly per dimension
    const int totSamples = 10000;

    // Take the dim'th root to get number of samples per dimension
    unsigned int pointsPerDim = std::pow(totSamples, 1.0/dim);

    // Clamp so that pointsPerDim >= 2
    if(pointsPerDim < 2) {
        pointsPerDim = 2;
    }

    return multi_linspace(dim, pointsPerDim);
}

std::vector<std::vector<double>> multi_linspace(int dim, unsigned int points_per_dim)
{
    auto start = std::vector<double>(dim);
    auto end = std::vector<double>(dim);
    auto numPoints = std::vector<unsigned int>(dim);


    for (int i = 0; i < dim; i++) {

        start.at(i) = -10.0;
        end.at(i) = 10.0;
        numPoints.at(i) = points_per_dim;
    }

    return multi_linspace(start, end, numPoints);
}

std::string pretty_print(const DenseVector &vec)
{
    std::string str("[");
    for (int i = 0; i < vec.rows(); i++) {
        str += to_string(vec(i));
        if(i + 1 < vec.rows()) {
            str += "; ";
        }
    }
    str += "]";

    return str;
}

} // namespace SPLINTER

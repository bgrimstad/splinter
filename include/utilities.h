/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_UTILITIES_H
#define SPLINTER_UTILITIES_H

#include <vector>
#include <cstdlib> // For std::abs, etc.
#include <definitions.h>

namespace SPLINTER
{

// Compare two numbers
template<typename T>
bool assert_near(T x, T y, double tolAbs = 1e-8, double tolRel = 1e-8)
{
    double dx = std::abs(x - y);
    double xAbs = 0.5*(std::abs(x) + std::abs(y));
    double err = std::max(tolAbs, tolRel*xAbs);
    return dx < err;
}

// Compare two vectors
inline bool compare_vectors(std::vector<double> x, std::vector<double> y, double tolAbs = 1e-8, double tolRel = 1e-8)
{
    if (x.size() != y.size())
        return false;

    for (unsigned int i = 0; i < x.size(); ++i)
        if (!assert_near(x.at(i), y.at(i), tolAbs, tolRel))
            return false;

    return true;
}

std::vector<double> eig_to_std_vec(const DenseVector &vec);

DenseVector std_to_eig_vec(const std::vector<double> &vec);

std::vector<std::vector<double>> eig_to_std_mat(const DenseMatrix &mat);

DenseMatrix std_to_eig_mat(const std::vector<std::vector<double>> &vec);

std::vector<double> linspace(double start, double stop, unsigned int num);

std::vector<double> extract_unique_sorted(const std::vector<double> &values);

std::vector<std::vector<double>> transpose_vec_vec(std::vector<std::vector<double>> x);

} // namespace SPLINTER

#endif // SPLINTER_UTILITIES_H

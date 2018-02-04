/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "utilities.h"

namespace SPLINTER
{

std::vector<double> eig_to_std_vec(const DenseVector &vec)
{
    std::vector<double> stdVec(vec.size());

    for(size_t i = 0; i < (size_t) stdVec.size(); ++i)
        stdVec.at(i) = vec(i);

    return stdVec;
}

DenseVector std_to_eig_vec(const std::vector<double> &vec)
{
    DenseVector eigVec(vec.size());
    eigVec.setZero();

    for (size_t i = 0; i < vec.size(); ++i)
        eigVec(i) = vec.at(i);

    return eigVec;
}

std::vector<std::vector<double>> eig_to_std_mat(const DenseMatrix &mat)
{
    std::vector<std::vector<double>> vec(mat.rows());

    for(size_t i = 0; i < (size_t) mat.rows(); ++i)
    {
        for(size_t j = 0; j < (size_t) mat.cols(); ++j)
        {
            vec.at(i).push_back(mat(i, j));
        }
    }

    return vec;
}

DenseMatrix std_to_eig_mat(const std::vector<std::vector<double>> &vec)
{
    size_t numRows = vec.size();
    size_t numCols = numRows > 0 ? vec.at(0).size() : 0;

    DenseMatrix mat(numRows, numCols);

    for(size_t i = 0; i < numRows; ++i)
    {
        for(size_t j = 0; j < numCols; ++j)
        {
            mat(i, j) = vec.at(i).at(j);
        }
    }

    return mat;
}

std::vector<double> linspace(double start, double stop, unsigned int num)
{
    std::vector<double> ret;
    double dx = 0;
    if (num > 1)
        dx = (stop - start)/(num-1);
    for (unsigned int i = 0; i < num; ++i)
        ret.push_back(start + i*dx);
    return ret;
}

std::vector<double> extract_unique_sorted(const std::vector<double> &values)
{
    // Sort and remove duplicates
    std::vector<double> unique(values);
    std::sort(unique.begin(), unique.end());
    std::vector<double>::iterator it = unique_copy(unique.begin(), unique.end(), unique.begin());
    unique.resize(distance(unique.begin(),it));
    return unique;
}

std::vector<std::vector<double>> transpose_vec_vec(std::vector<std::vector<double>> x)
{
    std::vector<std::vector<double>> xt;
    if (x.size() == 0) {
        return xt;
    }
    size_t inner_dim = x.at(0).size();

    for (unsigned int i = 0; i < inner_dim; ++i) {
        std::vector<double> inner;
        for (unsigned int j = 0; j < x.size(); ++j) {
            auto xji = x.at(j).at(i);
            inner.push_back(xji);
        }
        xt.push_back(inner);
    }
    return xt;
}

} // namespace SPLINTER

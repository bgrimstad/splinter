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

std::vector<double> denseVectorToVector(const DenseVector &denseVec)
{
    std::vector<double> vec(denseVec.size());

    for(size_t i = 0; i < denseVec.size(); ++i)
    {
        vec.at(i) = denseVec(i);
    }

    return vec;
}

DenseVector vectorToDenseVector(const std::vector<double> &vec)
{
    DenseVector denseVec(vec.size());

    for(size_t i = 0; i < vec.size(); ++i)
    {
        denseVec(i) = vec.at(i);
    }

    return denseVec;
}

std::vector<std::vector<double>> denseMatrixToVectorVector(const DenseMatrix &mat)
{
    std::vector<std::vector<double>> vec(mat.rows());

    for(size_t i = 0; i < mat.rows(); ++i)
    {
        vec.push_back(std::vector<double>());

        for(size_t j = 0; j < mat.cols(); ++j)
        {
            vec.at(i).push_back(mat(i, j));
        }
    }

    return vec;
}

DenseMatrix vectorVectorToDenseMatrix(const std::vector<std::vector<double>> &vec)
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

} // namespace SPLINTER
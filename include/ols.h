/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_OLS_H
#define SPLINTER_OLS_H

#include <definitions.h>
#include <linearfunction.h>
#include <datatable.h>
#include <linearsolvers.h>
#include <utilities.h>

namespace SPLINTER
{

/**
  * Ordinary least-square (OLS)
  * TODO: Select solver based on matrix type
  */
template <class Vec, class Mat>
DenseVector computeCoefficients(const LinearFunction<Vec, Mat> &func, const DataTable &data)
{
    // Left hand side
    Mat X = computeDesignMatrix(func, data);

    // Right-hand side
    auto y = vectorToDenseVector(data.getVectorY());

    // Coefficients
    DenseVector c;

    // Solve for coefficients
    DenseQR<> s;
    if (!s.solve(X, y, c))
        throw Exception("computeCoefficients: Failed to solve for coefficients.");

    return c;
};

/**
 * Computes design matrix by evaluating basis functions at each sample point
 */
template<class Vec, class Mat>
Mat computeDesignMatrix(const LinearFunction<Vec, Mat> &func, const DataTable &data);

} // namespace SPLINTER

#endif // SPLINTER_OLS_H

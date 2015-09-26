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
#include <sample.h>

namespace SPLINTER
{

/**
  * Ordinary least-square (OLS)
  */
DenseMatrix computeCoefficients(const LinearFunction &func, const Sample &samples);

// TODO: implement OLS with regularization term (lambda/numSample)*coefficients^T*coefficients,
// where lambda >= 0 and default is 0.1/N?
//DenseVector computeCoefficientsRegularized(const LinearFunction &func, const Sample &sample, double lambda = 0.1);

/**
 * Computes design matrix by evaluating basis functions
 * at each sample point
 */
DenseMatrix computeDesignMatrix(const LinearFunction &func, const Sample &sample);

SparseMatrix computeDesignMatrixSparse(const LinearFunction &func, const Sample &sample);

} // namespace SPLINTER

#endif // SPLINTER_OLS_H

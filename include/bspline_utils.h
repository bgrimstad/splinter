/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_UTILS_H
#define SPLINTER_BSPLINE_UTILS_H

#include <data_table.h>
#include <bspline.h>

namespace SPLINTER
{

// Matrix of basis functions evaluated at samples
SparseMatrix computeBasisFunctionMatrix(const BSpline &bspline, const DataTable &data);

DenseMatrix stackSamplePointValues(const DataTable &data);

// P-spline control point calculation
SparseMatrix computeSecondOrderFiniteDifferenceMatrix(const BSpline &bspline);

// Compute weights matrix from weight vector
SparseMatrix computeWeightMatrix(const std::vector<double> weights);

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_UTILS_H

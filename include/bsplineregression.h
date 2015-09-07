/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEREGRESSION_H
#define SPLINTER_BSPLINEREGRESSION_H

#include "datatable.h"
#include "bspline.h"

namespace SPLINTER
{

// Enum for different B-spline types
enum class BSplineType
{
    LINEAR,     // Linear basis functions in each variable
    QUADRATIC,  // Quadratic basis functions in each variable
    CUBIC,      // Cubic basis functions in each variable
    QUARTIC     // Quartic basis functions in each variable
};

BSpline buildBSpline(const DataTable &samples, std::vector<unsigned int> basisDegrees);
BSpline buildBSpline(const DataTable &samples, BSplineType type);

// Control point computations
DenseMatrix computeControlPoints(const DataTable &samples, const BSpline &bspline);
SparseMatrix computeBasisFunctionMatrix(const DataTable &samples, const BSpline &bspline);
DenseMatrix controlPointEquationRHS(const DataTable &samples);

// Computing knots
std::vector<std::vector<double> > computeKnotVectorsFromSamples(const DataTable &samples, std::vector<unsigned int> degrees);
std::vector<double> knotVectorMovingAverage(std::vector<double> &vec, unsigned int degree);
std::vector<double> knotVectorBuckets(std::vector<double> &vec, unsigned int degree);

} // namespace SPLINTER

#endif // SPLINTER_BSPLINEREGRESSION_H

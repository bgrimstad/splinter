/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_PSPLINE_H
#define SPLINTER_PSPLINE_H

#include "bsplineregression.h"

namespace SPLINTER
{

/*
 * The P-Spline is a smooting spline which relaxes the interpolation constraints on the control points to allow smoother spline curves.
 * It minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing).
 * It inherits all properties of the B-spline - the only difference lies in the calculation of the control points.
 */
BSpline computePSpline(const DataTable &samples, double lambda = 0.03);

// P-spline control point calculation
DenseMatrix computeControlPointsPSpline(const DataTable &samples, const BSpline &bspline, double lambda);
SparseMatrix getSecondOrderFiniteDifferenceMatrix(const BSpline &bspline);

} // namespace SPLINTER

#endif // SPLINTER_PSPLINE_H

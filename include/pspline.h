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
class SPLINTER_API PSpline : public BSplineRegression
{
public:
    PSpline(const char *fileName);
    PSpline(const std::string fileName);
    PSpline(const DataTable &samples);
    PSpline(const DataTable &samples, double lambda);

    double getLambda() { return lambda; }

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

protected:
    PSpline();

    // Smoothing parameter (usually set to a small number; default 0.03)
    double lambda;

    // P-spline control point calculation
    DenseMatrix computeControlPoints(const DataTable &samples) override;
    SparseMatrix getSecondOrderFiniteDifferenceMatrix();

private:
    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const PSpline &lhs, const PSpline &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_PSPLINE_H

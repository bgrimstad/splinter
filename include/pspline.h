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

#include "bspline.h"

namespace SPLINTER
{

/*
 * The P-Spline is a smooting spline which relaxes the interpolation constraints on the control points to allow smoother spline curves.
 * It minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing).
 * It inherits all properties of the B-spline - the only difference lies in the calculation of the control points.
 */
class API PSpline : public BSpline
{
public:
    PSpline();

    PSpline(const char *fileName);
    PSpline(const std::string fileName);
    PSpline(const DataTable &samples);
    PSpline(const DataTable &samples, double lambda);

    double getLambda() const { return lambda; }

    void save(const std::string fileName) const override;
    void _deserialize(StreamType::const_iterator &it, StreamType::const_iterator end);

protected:

    // Smoothing parameter (usually set to a small number; default 0.03)
    double lambda;

    // P-spline control point calculation
    void computeControlPoints(const DataTable &samples) override;
    void getSecondOrderFiniteDifferenceMatrix(SparseMatrix &D);

private:
    void load(const std::string fileName) override;
};

} // namespace SPLINTER

#endif // SPLINTER_PSPLINE_H

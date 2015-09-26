/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_PSPLINEAPPROXIMANT_H
#define SPLINTER_PSPLINEAPPROXIMANT_H

#include "bsplineapproximant.h"

namespace SPLINTER
{

/*
 * The P-Spline is a smooting spline which relaxes the interpolation constraints on the control points to allow smoother spline curves.
 * It minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing).
 * It inherits all properties of the B-spline - the only difference lies in the calculation of the control points.
 */
class SPLINTER_API PSplineApproximant : public BSplineApproximant
{
public:
    PSplineApproximant(const char *fileName);
    PSplineApproximant(const std::string fileName);
    PSplineApproximant(const Sample &samples, std::vector<unsigned int> basisDegrees, double lambda);
    PSplineApproximant(const Sample &samples, BSplineType type, double lambda);
    PSplineApproximant(const Sample &samples, double lambda);
    PSplineApproximant(const Sample &samples) : PSplineApproximant(samples, 0.03) {}

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

private:
    PSplineApproximant();

    double lambda; // Smoothing parameter. Requirement: lambda >= 0

    // P-spline control point calculation
    DenseMatrix computeCoefficients(const Sample &samples) const override;
    SparseMatrix getSecondOrderFiniteDifferenceMatrix() const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const PSplineApproximant &lhs, const PSplineApproximant &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_PSPLINEAPPROXIMANT_H

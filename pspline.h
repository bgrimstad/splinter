/*
This file is part of the Multivariate Splines library.
Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#ifndef PSPLINE_H
#define PSPLINE_H

#include "bspline.h"

namespace MultivariateSplines
{

/*
 * The P-Spline is a smooting spline which relaxes the interpolation constraints on the control points to allow smoother spline curves.
 * It minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing).
 * It inherits all properties of the B-spline - the only difference lies in the calculation of the control points.
 */
class PSpline : public BSpline
{
public:
    PSpline(DataTable &samples);
    PSpline(DataTable &samples, double lambda);

private:

    // Smoothing parameter (usually set to a small number; default 0.03)
    double lambda;

    // P-spline control point calculation
    void computeControlPoints(const DataTable &samples, std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y);
    void controlPointEquationLHS(const DataTable &samples, std::vector< std::vector<double> > &X, SparseMatrix &L, SparseMatrix &B, SparseMatrix &W, double lambda);
    void controlPointEquationRHS(std::vector< std::vector<double> > &Y, DenseMatrix &R, SparseMatrix &B, SparseMatrix &W);
    void getSecondOrderFiniteDifferenceMatrix(SparseMatrix &D);
    void getBasisFunctionMatrix(std::vector< std::vector<double> > &X, SparseMatrix &B);

};

} // namespace MultivariateSplines

#endif // PSPLINE_H

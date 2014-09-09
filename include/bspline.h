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


#ifndef BSPLINE_H
#define BSPLINE_H

#include "datatable.h"
#include "generaldefinitions.h"
#include "spline.h"
#include "basis.h"

namespace MultivariateSplines
{

// Enum for different B-spline types
enum class BSplineType
{
    LINEAR,             // Piecewise linear interpolation. Interpolates all points.
    //CUBIC_HERMITE,    // Cubic spline with Hermite end conditions. Interpolates all points. Not implemented.
    //CUBIC_NATURAL,    // Cubic spline with Natural end conditions. Interpolates all points. Ensures second derivative of B-spline is zero at end points. Not implemented.
    CUBIC_FREE          // Cubic spline with Free end conditions. Interpolates all points. Ensures p'th derivative continuous at x(2) and x(n-1). p+1-regular knot sequence with two deleted knots.
    //CUBIC_PERIODIC,   // Cubic spline with Periodic end conditions. Not implemented.
};

/*
 * Class that implements the multivariate tensor product B-spline
 */
class BSpline : public Spline
{
public:

    // Construct B-spline from knot sequences, control coefficients (assumed vectorized), and basis degrees
    //Bspline(std::vector<double> coefficients, std::vector<double> knotSequence, int basisDegrees);
    //Bspline(std::vector<double> coefficients, std::vector< std::vector<double> > knotSequences, std::vector<int> basisDegrees);
    BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotSequences, std::vector<int> basisDegrees);

    // Construct B-spline that interpolates the samples in DataTable
    BSpline(DataTable &samples, BSplineType type);

    virtual BSpline* clone() const { return new BSpline(*this); }

    void init();

    // Evaluation of B-spline
    double eval(DenseVector &x) const;
    DenseMatrix evalJacobian(DenseVector &x) const;
    DenseMatrix evalHessian(DenseVector &x) const;

    // Getters
    unsigned int getNumVariables() const { return numVariables; }
    unsigned int getNumControlPoints() const { return coefficients.cols(); }

    std::vector< std::vector<double> > getKnotVectors() const;

    std::vector<double> getDomainUpperBound() const;
    std::vector<double> getDomainLowerBound() const;

    // Control point related
    void setControlPoints(DenseMatrix &controlPoints);
    DenseMatrix getControlPoints() const;
    bool checkControlPoints() const;

    // B-spline operations
    bool reduceDomain(std::vector<double> lb, std::vector<double> ub, bool regularKnotsequences = true, bool refineKnotsequences = true);

    bool insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1); // TODO: move back to private

protected:

    BSpline() {}

    Basis basis;
    DenseMatrix knotaverages; // One row per input
    DenseMatrix coefficients; // One row per output

    unsigned int numVariables; // Dimension of x

    // Control point computations
    void computeKnotAverages();
    virtual void computeControlPoints(const DataTable &samples);
    void computeBasisFunctionMatrix(const DataTable &samples, SparseMatrix &A) const;
    void controlPointEquationRHS(const DataTable &samples, DenseMatrix &Bx, DenseMatrix &By) const;

private:

    // Domain reduction
    bool regularSequences(std::vector<double> &lb, std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub);

    // Knot insertion and refinement
    bool refineKnotSequences(); // All knots in one shabang

    // Helper functions
    bool valueInsideDomain(DenseVector x) const;

};

} // namespace MultivariateSplines

#endif // BSPLINE_H

/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef MS_BSPLINE_H
#define MS_BSPLINE_H

#include "datatable.h"
#include "generaldefinitions.h"
#include "spline.h"
#include "bsplinebasis.h"

namespace MultivariateSplines
{

// Enum for different B-spline types
enum class BSplineType
{
    LINEAR,             // Piecewise linear interpolation. Interpolates all points.
    QUADRATIC_FREE,     // Quadratic spline with free end conditions.
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

    // Construct B-spline from knot vectors, control coefficients (assumed vectorized), and basis degrees
    //Bspline(std::vector<double> coefficients, std::vector<double> knotVectors, unsigned int basisDegrees);
    //Bspline(std::vector<double> coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees);
    BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees);

    // Construct B-spline that interpolates the samples in DataTable
    //BSpline(DataTable &samples, unsigned int basisDegree);
    //BSpline(DataTable &samples, std::vector<unsigned int> basisDegrees);
    BSpline(const DataTable &samples, BSplineType type);

    virtual BSpline* clone() const { return new BSpline(*this); }

    void init();

    // Evaluation of B-spline
    double eval(DenseVector x) const;
    DenseMatrix evalJacobian(DenseVector x) const;
    DenseMatrix evalHessian(DenseVector x) const;

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
    bool reduceDomain(std::vector<double> lb, std::vector<double> ub, bool doRegularizeKnotVectors = true, bool doRefineKnotVectors = false);

    bool insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1); // TODO: move back to private

protected:

    BSpline() {}

    BSplineBasis basis;
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
    bool regularizeKnotVectors(std::vector<double> &lb, std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub);

    // Knot insertion and refinement
    bool refineKnotVectors(); // All knots in one shabang

    // Helper functions
    bool pointInDomain(DenseVector x) const;

};

} // namespace MultivariateSplines

#endif // MS_BSPLINE_H

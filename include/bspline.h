/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef SPLINTER_BSPLINE_H
#define SPLINTER_BSPLINE_H

#include "datatable.h"
#include "generaldefinitions.h"
#include "spline.h"
#include "bsplinebasis.h"

namespace SPLINTER
{

// Enum for different B-spline types
enum class BSplineType
{
    LINEAR,             // Piecewise linear interpolation. Interpolates all points.
    QUADRATIC,          // Quadratic spline with automatic selection of knots.
    QUADRATIC_FREE,     // Quadratic spline with free end conditions.
    CUBIC,              // Cubic spline with automatic selection of knots.
    //CUBIC_HERMITE,    // Cubic spline with Hermite end conditions. Interpolates all points. Not implemented.
    //CUBIC_NATURAL,    // Cubic spline with Natural end conditions. Interpolates all points. Ensures second derivative of B-spline is zero at end points. Not implemented.
    CUBIC_FREE          // Cubic spline with Free end conditions. Interpolates all points. Ensures p'th derivative continuous at x(2) and x(n-1).
    //CUBIC_PERIODIC,   // Cubic spline with Periodic end conditions. Not implemented.
};

/**
 * Class that implements the multivariate tensor product B-spline
 */
class BSpline : public Approximant
{
public:

    /**
     * Construct B-spline from knot vectors, control coefficients (assumed vectorized), and basis degrees
     */
    BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees);

    /**
     * Construct B-spline that interpolates the samples in DataTable
     */
    BSpline(const DataTable &samples, BSplineType type);

    /**
     * Construct B-spline from file
     */
    BSpline(const std::string fileName);

    virtual BSpline* clone() const { return new BSpline(*this); }

    void init();

    // Evaluation of B-spline
    double eval(DenseVector x) const;
    DenseMatrix evalJacobian(DenseVector x) const;
    DenseMatrix evalHessian(DenseVector x) const;

    // Getters
    unsigned int getNumVariables() const { return numVariables; }
    unsigned int getNumControlPoints() const { return coefficients.cols(); }
    std::vector<unsigned int> getNumBasisFunctions() const;
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<unsigned int> getBasisDegrees() const;
    std::vector<double> getDomainUpperBound() const;
    std::vector<double> getDomainLowerBound() const;

    // Control point related
    void setControlPoints(DenseMatrix &controlPoints);
    DenseMatrix getControlPoints() const;
    void checkControlPoints() const;

    // B-spline operations
    void reduceDomain(std::vector<double> lb, std::vector<double> ub, bool doRegularizeKnotVectors = true);

    // Perform global knot refinement
    void globalKnotRefinement(); // All knots in one shabang

    // Perform a local knot refinement at x
    void localKnotRefinement(DenseVector x);

    /**
     * Decompose B-spline to Bezier form
     */
    void decomposeToBezierForm();

    void insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1); // TODO: move back to private after testing

    void save(const std::string fileName) const override;

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
    void regularizeKnotVectors(std::vector<double> &lb, std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub);

    // Helper functions
    bool pointInDomain(DenseVector x) const;

    void load(const std::string fileName) override;
    void loadBasis(std::vector<std::vector<double>> knotVectors, std::vector<unsigned int> basisDegrees);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_H

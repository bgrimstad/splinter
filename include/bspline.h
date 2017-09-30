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

#include "function.h"
#include "bspline_basis.h"
#include "json_parser.h"


namespace SPLINTER
{

/**
 * Class that implements the multivariate tensor product B-spline
 */
class SPLINTER_API BSpline : public Function, public Saveable
{
public:
    /**
     * Builder class for construction by regression
     * Implemented in bspline_builder.*
     */
    class Builder;
    enum class Smoothing;

    /**
     * Construct B-spline from knot vectors, control points, and basis degrees
     */
    BSpline(unsigned int dimX,
            unsigned int dimY,
            const std::vector<std::vector<double>> &knotVectors,
            const std::vector<unsigned int> &degrees);

    BSpline(const std::vector<std::vector<double>> &controlPoints,
            const std::vector<std::vector<double>> &knotVectors,
            const std::vector<unsigned int> &degrees);

    /**
     * Construct B-spline from file
     */
    BSpline(const char *fileName);
    BSpline(const std::string &fileName);

    virtual BSpline* clone() const { return new BSpline(*this); }

    /**
     * Evaluation of B-spline
     */

    // Avoid name hiding
    using Function::eval;
    using Function::evalJacobian;

    /**
     * Returns the (dimY) function values at x
     */
    std::vector<double> eval(const std::vector<double> &x) const override;
    DenseVector eval(const DenseVector &x) const override;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    DenseMatrix evalJacobian(const DenseVector &x) const override;

    /**
     * Returns the (dimY x dimX x dimX) Hessian evaluated at x
     */
    std::vector<std::vector<std::vector<double>>> evalHessian(const std::vector<double> &x) const;

    // Evaluation of B-spline basis functions
    SparseVector evalBasis(const DenseVector &x) const;
    SparseMatrix evalBasisJacobian(const DenseVector &x) const;

    /**
     * Getters
     */
    DenseMatrix getControlPoints() const
    {
        return controlPoints;
    }

    unsigned int getNumControlPoints() const
    {
        return (unsigned int) controlPoints.rows();
    }

    std::vector<unsigned int> getNumBasisFunctionsPerVariable() const;

    unsigned int getNumBasisFunctions() const
    {
        return basis.getNumBasisFunctions();
    }

//    DenseMatrix getControlPoints() const;
    DenseMatrix getKnotAverages() const {
        return computeKnotAverages();
    };
    std::vector< std::vector<double>> getKnotVectors() const;
    std::vector<unsigned int> getBasisDegrees() const;
    std::vector<double> getDomainUpperBound() const;
    std::vector<double> getDomainLowerBound() const;

    /**
     * Setters
     */
    void setControlPoints(const DenseMatrix &newControlPoints);
    void checkControlPoints() const;

    // Linear transformation of control points (B-spline has affine invariance)
    void updateControlPoints(const SparseMatrix &A);

    // Reduce support of B-spline
    void reduceSupport(const std::vector<double> &lb, const std::vector<double> &ub, bool doRegularizeKnotVectors = true);

    // Perform global knot refinement
    void globalKnotRefinement(); // All knots in one shabang

    // Perform a local knot refinement at x
    void localKnotRefinement(const DenseVector &x);

    // Decompose B-spline to Bezier form
    void decomposeToBezierForm();

    // Insert a knot until desired knot multiplicity is obtained
    void insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1);

    /**
     * Save and load
     */
    void save(const std::string &fileName) const override;

    void save_to_json(const std::string &filename) const;

    static BSpline load_from_json(const std::string &filename);

    std::string getDescription() const override;

protected:
    BSpline();

    BSplineBasis basis;

    /*
     * B-spline control points in R^(m x n),
     * where m = numBasisFunctions and n = dim_y.
     * Each row is a control point.
     */
    DenseMatrix controlPoints;

    // Control point computations
    DenseMatrix computeKnotAverages() const;

private:
    // Domain reduction
    void regularizeKnotVectors(const std::vector<double> &lb, const std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(const std::vector<double> &lb, const std::vector<double> &ub);

    // Helper functions
    bool isSupported(const DenseVector &x) const;

    void load(const std::string &fileName) override;

    friend class Serializer;
    friend bool operator==(const BSpline &lhs, const BSpline &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_H

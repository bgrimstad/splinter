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
#include "function.h"
#include "bsplinebasis.h"

namespace SPLINTER
{

/**
 * Class that implements the multivariate tensor product B-spline
 */
class SPLINTER_API BSpline : public Function
{
public:
    BSpline(unsigned int numVariables);

    /**
     * Construct B-spline from knot vectors, control coefficients (assumed vectorized), and basis degrees
     */
    BSpline(std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees); // All coefficients set to 1
    BSpline(std::vector<double> coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees);
    BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees);

    /**
     * Construct B-spline from file
     */
    BSpline(const char *fileName);
    BSpline(const std::string fileName);

    virtual BSpline* clone() const { return new BSpline(*this); }

    // Evaluation of B-spline
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override;

    // Evaluation of B-spline basis functions
    SparseVector evalBasisFunctions(DenseVector x) const;
    SparseMatrix evalBasisFunctionsJacobian(DenseVector x) const;

    // Getters
    unsigned int getNumControlPoints() const { return coefficients.cols(); }
    std::vector<unsigned int> getNumBasisFunctions() const;
    unsigned int getNumBasisFunctionsTotal() const
    {
        return basis.getNumBasisFunctions();
    }
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<unsigned int> getBasisDegrees() const;
    std::vector<double> getDomainUpperBound() const;
    std::vector<double> getDomainLowerBound() const;

    // Control point related
    DenseMatrix getControlPoints() const;
    void updateControlPoints(const DenseMatrix &A);
    void setCoefficients(const DenseMatrix &coefficients);
    void setControlPoints(const DenseMatrix &controlPoints);
    void checkControlPoints() const;

    // B-spline operations
    void reduceDomain(std::vector<double> lb, std::vector<double> ub, bool doRegularizeKnotVectors = true);

    // Perform global knot refinement
    void globalKnotRefinement(); // All knots in one shabang

    // Perform a local knot refinement at x
    void localKnotRefinement(DenseVector x);

    // Decompose B-spline to Bezier form
    void decomposeToBezierForm();

    void insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1); // TODO: move back to private after testing

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

protected:
    BSpline();

    BSplineBasis basis;
    DenseMatrix knotaverages; // One row per input
    DenseMatrix coefficients; // One row per output

    // Control point computations
    void computeKnotAverages();

private:
    // Domain reduction
    void regularizeKnotVectors(std::vector<double> &lb, std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub);

    // Helper functions
    bool pointInDomain(DenseVector x) const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const BSpline &lhs, const BSpline &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_H

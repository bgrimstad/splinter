/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_BASIS_H
#define SPLINTER_BSPLINE_BASIS_H

#include "definitions.h"
#include "bspline_basis_1d.h"

namespace SPLINTER
{

class BSplineBasis
{
public:
    BSplineBasis();
    BSplineBasis(const std::vector<std::vector<double>> &knotVectors, std::vector<unsigned int> basisDegrees);

    // Evaluation
    SparseVector eval(const DenseVector &x) const;
    DenseMatrix evalBasisJacobianOld(const DenseVector &x) const; // Deprecated
    SparseMatrix evalBasisJacobian(const DenseVector &x) const;
    SparseMatrix evalBasisJacobian2(const DenseVector &x) const; // A bit slower than evaBasisJacobianOld()
    SparseMatrix evalBasisHessian(const DenseVector &x) const;

    // Knot vector manipulation
    SparseMatrix refineKnots();
    SparseMatrix refineKnotsLocally(const DenseVector &x);
    SparseMatrix decomposeToBezierForm();
    SparseMatrix insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1);

    // Getters
    BSplineBasis1D getSingleBasis(unsigned int dim);
    std::vector<std::vector<double>> getKnotVectors() const;
    std::vector<double> getKnotVector(int dim) const;

    std::vector<unsigned int> getBasisDegrees() const;
    unsigned int getBasisDegree(unsigned int dim) const;
    unsigned int getNumBasisFunctions() const;
    unsigned int getNumBasisFunctions(unsigned int dim) const;
    std::vector<unsigned int> getNumBasisFunctionsTarget() const;

    unsigned int getKnotMultiplicity(unsigned int dim, double tau) const;

    /*
     * Returns the maximum number of supported basis functions at any point in the B-spline domain
     */
    unsigned int numSupported() const;

    bool insideSupport(const DenseVector &x) const;
    std::vector<double> getSupportLowerBound() const;
    std::vector<double> getSupportUpperBound() const;

    // Support related
    SparseMatrix reduceSupport(const std::vector<double>& lb, const std::vector<double>& ub);

private:
    std::vector<BSplineBasis1D> bases;
    unsigned int numVariables;

    friend class Serializer;
    friend bool operator==(const BSplineBasis &lhs, const BSplineBasis &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BASIS_H

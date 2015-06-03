/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEBASIS_H
#define SPLINTER_BSPLINEBASIS_H

#include "generaldefinitions.h"
#include "bsplinebasis1d.h"

namespace SPLINTER
{

class BSplineBasis
{
public:
    BSplineBasis();
    BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees);
    BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees, KnotVectorType knotVectorType);

    void setUnivariateBases(std::vector< std::vector<double> > &X, std::vector<unsigned int> &basisDegrees, KnotVectorType knotVectorType);

    // Evaluation
    SparseVector eval(const DenseVector &x) const;
    SparseMatrix evalBasisJacobian(DenseVector &x) const;
    DenseMatrix evalBasisJacobianOld(DenseVector &x) const; // Depricated
    DenseMatrix evalBasisJacobianFast(DenseVector &x) const; // Depricated
    SparseMatrix evalBasisHessian(DenseVector &x) const;

    // Knot vector manipulation
    SparseMatrix refineKnots();
    SparseMatrix refineKnotsLocally(DenseVector x);
    SparseMatrix decomposeToBezierForm();
    SparseMatrix insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1);

    // Getters
    BSplineBasis1D getSingleBasis(int dim);
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<double> getKnotVector(int dim) const;

    std::vector<unsigned int> getBasisDegrees() const;
    unsigned int getBasisDegree(unsigned int dim) const;
    unsigned int getNumBasisFunctions() const;
    unsigned int getNumBasisFunctions(unsigned int dim) const;
    std::vector<unsigned int> getNumBasisFunctionsTarget() const;

    double getKnotValue(int dim, int index) const;
    unsigned int getKnotMultiplicity(unsigned int dim, double tau) const;
    unsigned int getLargestKnotInterval(unsigned int dim) const;

    int supportedPrInterval() const;

    bool insideSupport(DenseVector &x) const;
    std::vector<double> getSupportLowerBound() const;
    std::vector<double> getSupportUpperBound() const;

    // Support related
    bool reduceSupport(std::vector<double>& lb, std::vector<double>& ub, SparseMatrix &A);

private:
    std::vector<BSplineBasis1D> bases;
    unsigned int numVariables;
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINEBASIS_H

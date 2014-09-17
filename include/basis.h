/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef BASIS_H
#define BASIS_H

#include "generaldefinitions.h"
#include "basis1d.h"

namespace MultivariateSplines
{

class Basis
{
public:
    Basis();
    Basis(std::vector< std::vector<double> > &X, std::vector<int> basisDegrees);
    Basis(std::vector< std::vector<double> > &X, std::vector<int> basisDegrees, KnotSequenceType knotSequenceType);

    void setUnivariateBases(std::vector< std::vector<double> > &X, std::vector<int> &basisDegrees, KnotSequenceType knotSequenceType);

    // Evaluation
    SparseVector eval(const DenseVector &x) const;
    SparseMatrix evalBasisJacobian(DenseVector &x) const;
    DenseMatrix evalBasisJacobianOld(DenseVector &x) const; // Depricated
    SparseMatrix evalBasisHessian(DenseVector &x) const;

    // Knot insertion
    bool refineKnots(SparseMatrix &A);
    bool insertKnots(SparseMatrix &A, double tau, unsigned int dim, unsigned int multiplicity = 1);
    //bool insertKnots(SparseMatrix &A, std::vector<std::tuple<double,int,int>> tau, unsigned int dim, unsigned int multiplicity = 1);

    // Getters
    Basis1D getSingleBasis(int dim);
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<double> getKnotVector(int dim) const;

    int getBasisDegree(int dim) const;
    int numBasisFunctions() const;
    int numBasisFunctions(int dim) const;

    double getKnotValue(int dim, int index) const;
    int getKnotMultiplicity(const int& dim, const double &tau) const;
    int getLargestKnotInterval(int dim) const;

    std::vector<int> getTensorIndexDimension() const;
    std::vector<int> getTensorIndexDimensionTarget() const;
    int supportedPrInterval() const;

    bool valueInsideSupport(DenseVector &x) const;
    std::vector<double> getSupportLowerBound() const;
    std::vector<double> getSupportUpperBound() const;

    // Support related
    bool reduceSupport(std::vector<double>& lb, std::vector<double>& ub, SparseMatrix &A);

private:
    std::vector<Basis1D> bases;
    unsigned int numVariables;
};

} // namespace MultivariateSplines

#endif // BASIS_H

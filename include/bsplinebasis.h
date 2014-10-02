/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef MS_BSPLINEBASIS_H
#define MS_BSPLINEBASIS_H

#include "generaldefinitions.h"
#include "bsplinebasis1d.h"

namespace MultivariateSplines
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
    SparseMatrix evalBasisHessian(DenseVector &x) const;

    // Knot insertion
    bool refineKnots(SparseMatrix &A);
    bool insertKnots(SparseMatrix &A, double tau, unsigned int dim, unsigned int multiplicity = 1);
    //bool insertKnots(SparseMatrix &A, std::vector<std::tuple<double,int,int>> tau, unsigned int dim, unsigned int multiplicity = 1);

    // Getters
    BSplineBasis1D getSingleBasis(int dim);
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<double> getKnotVector(int dim) const;

    unsigned int getBasisDegree(unsigned int dim) const;
    unsigned int numBasisFunctions() const;
    unsigned int numBasisFunctions(unsigned int dim) const;

    double getKnotValue(int dim, int index) const;
    unsigned int getKnotMultiplicity(unsigned int dim, double tau) const;
    unsigned int getLargestKnotInterval(unsigned int dim) const;

    std::vector<int> getTensorIndexDimension() const;
    std::vector<int> getTensorIndexDimensionTarget() const;
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

} // namespace MultivariateSplines

#endif // MS_BSPLINEBASIS_H

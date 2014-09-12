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


#ifndef MS_BSPLINEBASIS_H
#define MS_BSPLINEBASIS_H

#include "generaldefinitions.h"
#include "basis1d.h"

namespace MultivariateSplines
{

class Basis
{
public:
    Basis();
    Basis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees);
    Basis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees, KnotVectorType knotVectorType);

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
    Basis1D getSingleBasis(int dim);
    std::vector< std::vector<double> > getKnotVectors() const;
    std::vector<double> getKnotVector(int dim) const;

    unsigned int getBasisDegree(unsigned int dim) const;
    unsigned int numBasisFunctions() const;
    unsigned int numBasisFunctions(unsigned int dim) const;

    double getKnotValue(int dim, int index) const;
    unsigned int getKnotMultiplicity(const unsigned int& dim, const double &tau) const;
    unsigned int getLargestKnotInterval(unsigned int dim) const;

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

#endif // MS_BSPLINEBASIS_H

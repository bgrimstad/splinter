/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "include/bsplinebasis.h"
#include "include/mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"

namespace MultivariateSplines
{

BSplineBasis::BSplineBasis()
{
}

BSplineBasis::BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees)
    : BSplineBasis(X, basisDegrees, KnotVectorType::FREE)
{
}

BSplineBasis::BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees, KnotVectorType knotVectorType)
{
    assert(X.size() == basisDegrees.size());

    setUnivariateBases(X, basisDegrees, knotVectorType);
}

void BSplineBasis::setUnivariateBases(std::vector< std::vector<double> > &X, std::vector<unsigned int> &basisDegrees, KnotVectorType knotVectorType)
{
    bases.clear();
    numVariables = X.size();
    for(unsigned int i = 0; i < numVariables; i++)
    {
        bases.push_back(BSplineBasis1D(X.at(i), basisDegrees.at(i), knotVectorType));
    }
}

SparseVector BSplineBasis::eval(const DenseVector &x) const
{
    // Set correct starting size and values: one element of 1.
    SparseMatrix tv1(1,1);  //tv1.resize(1);
    tv1.reserve(1);
    tv1.insert(0,0) = 1;
    SparseMatrix tv2 = tv1;

    // Evaluate all basisfunctions in each dimension i and generate the tensor product of the function values
    for(int dim = 0; dim < x.size(); dim++)
    {
        SparseVector xi = bases.at(dim).evaluate(x(dim));

        // Avoid matrix copy
        if(dim % 2 == 0)
        {
            tv1 = kroneckerProduct(tv2, xi); // tv1 = tv1 x xi
        }
        else
        {
            tv2 = kroneckerProduct(tv1, xi); // tv2 = tv1 x xi
        }
    }

    // Return correct vector
    if(tv1.size() == numBasisFunctions())
    {
        return tv1;
    }
    return tv2;
}

// Old implementation of Jacobian
DenseMatrix BSplineBasis::evalBasisJacobianOld(DenseVector &x) const
{
    // Jacobian basis matrix
    DenseMatrix J; J.setZero(numBasisFunctions(), numVariables);

    // Calculate partial derivatives
    for(unsigned int i = 0; i < numVariables; i++)
    {
        // One column in basis jacobian
        DenseVector bi; bi.setOnes(1);
        for(unsigned int j = 0; j < numVariables; j++)
        {
            DenseVector temp = bi;
            DenseVector xi;
            if(j == i)
            {
                // Differentiated basis
                xi = bases.at(j).evaluateFirstDerivative(x(j));
            }
            else
            {
                // Normal basis
                xi = bases.at(j).evaluate(x(j));
            }

            bi = kroneckerProduct(temp, xi);
            //bi = kroneckerProduct(bi, xi).eval();
        }

        // Fill out column
        J.block(0,i,bi.rows(),1) = bi.block(0,0,bi.rows(),1);
    }

    return J;
}

SparseMatrix BSplineBasis::evalBasisJacobian(DenseVector &x) const
{
    // Jacobian basis matrix
    SparseMatrix J(numBasisFunctions(), numVariables);
    //J.setZero(numBasisFunctions(), numInputs);

    // Calculate partial derivatives
    for(unsigned int i = 0; i < numVariables; i++)
    {
        // One column in basis jacobian
        SparseMatrix Ji(1,1);
        Ji.insert(0,0) = 1;
        for(unsigned int j = 0; j < numVariables; j++)
        {
            SparseMatrix temp = Ji;
            SparseMatrix xi;
            if(j == i)
            {
                // Differentiated basis
                xi = bases.at(j).evaluateDerivative(x(j),1);
            }
            else
            {
                // Normal basis
                xi = bases.at(j).evaluate(x(j));
            }
            Ji = kroneckerProduct(temp, xi);
            //myKroneckerProduct(temp,xi,Ji);
            //Ji = kroneckerProduct(Ji, xi).eval();
        }

        // Fill out column
        for(int k = 0; k < Ji.outerSize(); ++k)
        for(SparseMatrix::InnerIterator it(Ji,k); it; ++it)
        {
            if(it.value() != 0)
                J.insert(it.row(),i) = it.value();
        }
        //J.block(0,i,Ji.rows(),1) = bi.block(0,0,Ji.rows(),1);
    }

    J.makeCompressed();

    return J;
}

SparseMatrix BSplineBasis::evalBasisHessian(DenseVector &x) const
{
    // Hessian basis matrix
    /* Hij = B1 x ... x DBi x ... x DBj x ... x Bn
     * (Hii = B1 x ... x DDBi x ... x Bn)
     * Where B are basis functions evaluated at x,
     * DB are the derivative of the basis functions,
     * and x is the kronecker product.
     * Hij is in R^(numBasisFunctions x 1)
     * so that basis hessian H is in R^(numBasisFunctions*numInputs x numInputs)
     * The real B-spline Hessian is calculated as (c^T x 1^(numInputs x 1))*H
     */
    SparseMatrix H(numBasisFunctions()*numVariables, numVariables);
    //H.setZero(numBasisFunctions()*numInputs, numInputs);

    // Calculate partial derivatives
    // Utilizing that Hessian is symmetric
    // Filling out lower left triangular
    for(unsigned int i = 0; i < numVariables; i++) // row
    {
        for(unsigned int j = 0; j <= i; j++) // col
        {
            // One column in basis jacobian
            SparseMatrix Hi(1,1);
            Hi.insert(0,0) = 1;

            for(unsigned int k = 0; k < numVariables; k++)
            {
                SparseMatrix temp = Hi;
                SparseMatrix Bk;
                if(i == j && k == i)
                {
                    // Diagonal element
                    Bk = bases.at(k).evaluateDerivative(x(k),2);
                }
                else if(k == i || k == j)
                {
                    Bk = bases.at(k).evaluateDerivative(x(k),1);
                }
                else
                {
                    Bk = bases.at(k).evaluate(x(k));
                }
                Hi = kroneckerProduct(temp, Bk);
            }

            // Fill out column
            for(int k = 0; k < Hi.outerSize(); ++k)
            for(SparseMatrix::InnerIterator it(Hi,k); it; ++it)
            {
                if(it.value() != 0)
                {
                    int row = i*numBasisFunctions()+it.row();
                    int col = j;
                    H.insert(row,col) = it.value();
                }
            }
        }
    }

    H.makeCompressed();

    return H;
}

bool BSplineBasis::insertKnots(SparseMatrix &A, double tau, unsigned int dim, unsigned int multiplicity)
{
    A.resize(1,1);
    A.insert(0,0) = 1;

    // Calculate multivariate knot insertion matrix
    for(unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        if(i == dim)
        {
            // Build knot insertion matrix
            if (!bases.at(i).insertKnots(Ai, tau, multiplicity))
                return false;
        }
        else
        {
            // No insertion - identity matrix
            int m = bases.at(i).numBasisFunctions();
            Ai.resize(m,m);
            Ai.setIdentity();
        }

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp,Ai,A);
    }

    A.makeCompressed();

    return true;
}

bool BSplineBasis::refineKnots(SparseMatrix &A)
{
    A.resize(1,1);
    A.insert(0,0) = 1;

    for(unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        if(!bases.at(i).refineKnots(Ai))
            return false;

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp,Ai,A);
    }

    A.makeCompressed();

    return true;
}

bool BSplineBasis::reduceSupport(std::vector<double>& lb, std::vector<double>& ub, SparseMatrix &A)
{
    assert(lb.size() == ub.size());
    assert(lb.size() == numVariables);

    A.resize(1,1);
    A.insert(0,0) = 1;

    for(unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        if(!bases.at(i).reduceSupport(lb.at(i), ub.at(i), Ai))
            return false;

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return true;
}

unsigned int BSplineBasis::getBasisDegree(unsigned int dim) const
{
    return bases.at(dim).getBasisDegree();
}

unsigned int BSplineBasis::numBasisFunctions(unsigned int dim) const
{
    return bases.at(dim).numBasisFunctions();
}

unsigned int BSplineBasis::numBasisFunctions() const
{
    unsigned int prod = 1;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        prod *= bases.at(dim).numBasisFunctions();
    }
    return prod;
}

BSplineBasis1D BSplineBasis::getSingleBasis(int dim)
{
    return bases.at(dim);
}

std::vector<double> BSplineBasis::getKnotVector(int dim) const
{
    return bases.at(dim).getKnotVector();
}

std::vector< std::vector<double> > BSplineBasis::getKnotVectors() const
{
    std::vector< std::vector<double> > knots;
    for(unsigned int i = 0; i < numVariables; i++)
        knots.push_back(bases.at(i).getKnotVector());
    return knots;
}

unsigned int BSplineBasis::getKnotMultiplicity(unsigned int dim, double tau) const
{
    return bases.at(dim).knotMultiplicity(tau);
}

double BSplineBasis::getKnotValue(int dim, int index) const
{
    return bases.at(dim).getKnotValue(index);
}

unsigned int BSplineBasis::getLargestKnotInterval(unsigned int dim) const
{
    return bases.at(dim).indexLongestInterval();
}

std::vector<int> BSplineBasis::getTensorIndexDimension() const
{
    std::vector<int> ret;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        ret.push_back(bases.at(dim).numBasisFunctions());
    }
    return ret;
}

std::vector<int> BSplineBasis::getTensorIndexDimensionTarget() const
{
    std::vector<int> ret;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        ret.push_back(bases.at(dim).numBasisFunctionsTarget() );
    }
    return ret;
}

int BSplineBasis::supportedPrInterval() const
{
    int ret = 1;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        ret *= (bases.at(dim).getBasisDegree() + 1);
    }
    return ret;
}

bool BSplineBasis::insideSupport(DenseVector &x) const
{
    if(x.size() != numVariables)
    {
        return false;
    }

    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        if(!bases.at(dim).insideSupport(x(dim)))
        {
            return false;
        }
    }
    return true;
}

std::vector<double> BSplineBasis::getSupportLowerBound() const
{
    std::vector<double> lb;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        std::vector<double> knots = bases.at(dim).getKnotVector();
        lb.push_back(knots.front());
    }
    return lb;
}

std::vector<double> BSplineBasis::getSupportUpperBound() const
{
    std::vector<double> ub;
    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        std::vector<double> knots = bases.at(dim).getKnotVector();
        ub.push_back(knots.back());
    }
    return ub;
}

} // namespace MultivariateSplines

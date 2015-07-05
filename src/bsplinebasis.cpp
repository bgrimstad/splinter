/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bsplinebasis.h"
#include "mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"

#include <iostream>

namespace SPLINTER
{

BSplineBasis::BSplineBasis()
{
}

BSplineBasis::BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees)
    : BSplineBasis(X, basisDegrees, false)
{
}

BSplineBasis::BSplineBasis(std::vector< std::vector<double> > &X, std::vector<unsigned int> basisDegrees, bool explicitKnots)
{
    if (X.size() != basisDegrees.size())
        throw Exception("BSplineBasis::BSplineBasis: Incompatible sizes. Number of knot vectors is not equal to size of degree vector.");

    setUnivariateBases(X, basisDegrees, explicitKnots);
}

BSplineBasis::BSplineBasis(std::vector<BSplineBasis1D> &bases, unsigned int numVariables)
    : bases(bases),
    numVariables(numVariables)
{
}

void BSplineBasis::setUnivariateBases(std::vector< std::vector<double> > &X, std::vector<unsigned int> &basisDegrees, bool explicitKnots)
{
    bases.clear();
    numVariables = X.size();
    for (unsigned int i = 0; i < numVariables; i++)
    {
        bases.push_back(BSplineBasis1D(X.at(i), basisDegrees.at(i), explicitKnots));

        // Adjust target number of basis functions used in e.g. refinement
        if (numVariables > 2)
        {
            // One extra knot is allowed
            bases.at(i).setNumBasisFunctionsTarget((basisDegrees.at(i)+1)+1); // Minimum degree+1
        }
    } 
}

SparseVector BSplineBasis::eval(const DenseVector &x) const
{
    // Evaluate basisfunctions for each variable i and compute the tensor product of the function values
    std::vector<SparseVector> basisFunctionValues;

    for (int var = 0; var < x.size(); var++)
        basisFunctionValues.push_back(bases.at(var).evaluate(x(var)));

    return foldlKroneckerProductVectors(basisFunctionValues);
}

// Old implementation of Jacobian
DenseMatrix BSplineBasis::evalBasisJacobianOld(DenseVector &x) const
{
    // Jacobian basis matrix
    DenseMatrix J; J.setZero(getNumBasisFunctions(), numVariables);

    // Calculate partial derivatives
    for (unsigned int i = 0; i < numVariables; i++)
    {
        // One column in basis jacobian
        DenseVector bi; bi.setOnes(1);
        for (unsigned int j = 0; j < numVariables; j++)
        {
            DenseVector temp = bi;
            DenseVector xi;
            if (j == i)
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
        }

        // Fill out column
        J.block(0,i,bi.rows(),1) = bi.block(0,0,bi.rows(),1);
    }

    return J;
}

SparseMatrix BSplineBasis::evalBasisJacobian(DenseVector &x) const
{
    // Jacobian basis matrix
    SparseMatrix J(getNumBasisFunctions(), numVariables);
    //J.setZero(numBasisFunctions(), numInputs);

    // Calculate partial derivatives
    for (unsigned int i = 0; i < numVariables; ++i)
    {
        // One column in basis jacobian
        SparseMatrix Ji(1,1);
        Ji.insert(0,0) = 1;
        for (unsigned int j = 0; j < numVariables; ++j)
        {
            SparseMatrix temp = Ji;
            SparseMatrix xi;
            if (j == i)
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
        for (int k = 0; k < Ji.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(Ji,k); it; ++it)
        {
            if (it.value() != 0)
                J.insert(it.row(),i) = it.value();
        }
        //J.block(0,i,Ji.rows(),1) = bi.block(0,0,Ji.rows(),1);
    }

    J.makeCompressed();

    return J;
}

SparseMatrix BSplineBasis::evalBasisJacobian2(DenseVector &x) const
{
    // Jacobian basis matrix
    SparseMatrix J(getNumBasisFunctions(), numVariables);

    // Evaluate B-spline basis functions before looping
    std::vector<SparseVector> funcValues(numVariables);
    std::vector<SparseVector> gradValues(numVariables);

    for (unsigned int i = 0; i < numVariables; ++i)
    {
        funcValues[i] = bases.at(i).evaluate(x(i));
        gradValues[i] = bases.at(i).evaluateFirstDerivative(x(i));
    }

    // Calculate partial derivatives
    for (unsigned int i = 0; i < numVariables; i++)
    {
        std::vector<SparseVector> values;

        for (unsigned int j = 0; j < numVariables; j++)
        {
            if (j == i)
                values.push_back(gradValues.at(j)); // Differentiated basis
            else
                values.push_back(funcValues.at(j)); // Normal basis
        }

        SparseVector Ji = foldlKroneckerProductVectors(values);

        // Fill out column
        for (SparseVector::InnerIterator it(Ji); it; ++it)
            J.insert(it.row(),i) = it.value();
    }

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
    SparseMatrix H(getNumBasisFunctions()*numVariables, numVariables);
    //H.setZero(numBasisFunctions()*numInputs, numInputs);

    // Calculate partial derivatives
    // Utilizing that Hessian is symmetric
    // Filling out lower left triangular
    for (unsigned int i = 0; i < numVariables; i++) // row
    {
        for (unsigned int j = 0; j <= i; j++) // col
        {
            // One column in basis jacobian
            SparseMatrix Hi(1,1);
            Hi.insert(0,0) = 1;

            for (unsigned int k = 0; k < numVariables; k++)
            {
                SparseMatrix temp = Hi;
                SparseMatrix Bk;
                if (i == j && k == i)
                {
                    // Diagonal element
                    Bk = bases.at(k).evaluateDerivative(x(k),2);
                }
                else if (k == i || k == j)
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
            for (int k = 0; k < Hi.outerSize(); ++k)
            for (SparseMatrix::InnerIterator it(Hi,k); it; ++it)
            {
                if (it.value() != 0)
                {
                    int row = i*getNumBasisFunctions()+it.row();
                    int col = j;
                    H.insert(row,col) = it.value();
                }
            }
        }
    }

    H.makeCompressed();

    return H;
}

SparseMatrix BSplineBasis::insertKnots(double tau, unsigned int dim, unsigned int multiplicity)
{
    SparseMatrix A(1,1);
//    A.resize(1,1);
    A.insert(0,0) = 1;

    // Calculate multivariate knot insertion matrix
    for (unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        if (i == dim)
        {
            // Build knot insertion matrix
            Ai = bases.at(i).insertKnots(tau, multiplicity);
        }
        else
        {
            // No insertion - identity matrix
            int m = bases.at(i).getNumBasisFunctions();
            Ai.resize(m,m);
            Ai.setIdentity();
        }

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::refineKnots()
{
    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).refineKnots();

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::refineKnotsLocally(DenseVector x)
{
    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).refineKnotsLocally(x(i));

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::decomposeToBezierForm()
{
    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai = bases.at(i).decomposeToBezierForm();

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return A;
}

SparseMatrix BSplineBasis::reduceSupport(std::vector<double>& lb, std::vector<double>& ub)
{
    assert(lb.size() == ub.size());
    assert(lb.size() == numVariables);

    SparseMatrix A(1,1);
    A.insert(0,0) = 1;

    for (unsigned int i = 0; i < numVariables; i++)
    {
        SparseMatrix temp = A;
        SparseMatrix Ai;

        Ai = bases.at(i).reduceSupport(lb.at(i), ub.at(i));

        //A = kroneckerProduct(temp, Ai);
        myKroneckerProduct(temp, Ai, A);
    }

    A.makeCompressed();

    return A;
}

SparseVector BSplineBasis::foldlKroneckerProductVectors(const std::vector<SparseVector> &vectors) const
{
    // Create two temp matrices
    SparseMatrix temp1(1,1);
    temp1.insert(0,0) = 1;
    SparseMatrix temp2 = temp1;

    // Fold from left
    int counter = 0;
    for (const auto &vec : vectors)
    {
        // Avoid copy
        if (counter % 2 == 0)
            temp1 = kroneckerProduct(temp2, vec);
        else
            temp2 = kroneckerProduct(temp1, vec);

        ++counter;
    }

    // Return correct product
    if (counter % 2 == 0)
        return temp2;
    return temp1;
}

SparseMatrix BSplineBasis::foldlKroneckerProductMatrices(const std::vector<SparseMatrix> &matrices) const
{
    // Create two temp matrices
    SparseMatrix temp1(1,1);
    temp1.insert(0,0) = 1;
    SparseMatrix temp2 = temp1;

    // Fold from left
    int counter = 0;
    for (const auto &mat : matrices)
    {
        // Avoid copy
        if (counter % 2 == 0)
            temp1 = kroneckerProduct(temp2, mat);
        else
            temp2 = kroneckerProduct(temp1, mat);

        ++counter;
    }

    // Return correct product
    if (counter % 2 == 0)
        return temp2;
    return temp1;
}

std::vector<unsigned int> BSplineBasis::getBasisDegrees() const
{
    std::vector<unsigned int> degrees;
    for (const auto& basis : bases)
        degrees.push_back(basis.getBasisDegree());
    return degrees;
}

unsigned int BSplineBasis::getBasisDegree(unsigned int dim) const
{
    return bases.at(dim).getBasisDegree();
}

unsigned int BSplineBasis::getNumBasisFunctions(unsigned int dim) const
{
    return bases.at(dim).getNumBasisFunctions();
}

unsigned int BSplineBasis::getNumBasisFunctions() const
{
    unsigned int prod = 1;
    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        prod *= bases.at(dim).getNumBasisFunctions();
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
    for (unsigned int i = 0; i < numVariables; i++)
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

std::vector<unsigned int> BSplineBasis::getNumBasisFunctionsTarget() const
{
    std::vector<unsigned int> ret;
    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        ret.push_back(bases.at(dim).getNumBasisFunctionsTarget() );
    }
    return ret;
}

int BSplineBasis::supportedPrInterval() const
{
    int ret = 1;
    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        ret *= (bases.at(dim).getBasisDegree() + 1);
    }
    return ret;
}

bool BSplineBasis::insideSupport(DenseVector &x) const
{
    if (x.size() != numVariables)
    {
        return false;
    }

    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        if (!bases.at(dim).insideSupport(x(dim)))
        {
            return false;
        }
    }
    return true;
}

std::vector<double> BSplineBasis::getSupportLowerBound() const
{
    std::vector<double> lb;
    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        std::vector<double> knots = bases.at(dim).getKnotVector();
        lb.push_back(knots.front());
    }
    return lb;
}

std::vector<double> BSplineBasis::getSupportUpperBound() const
{
    std::vector<double> ub;
    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        std::vector<double> knots = bases.at(dim).getKnotVector();
        ub.push_back(knots.back());
    }
    return ub;
}

std::vector<BSplineBasis1D> BSplineBasis::getBases() const
{
    return bases;
}

unsigned int BSplineBasis::getNumVariables() const
{
    return numVariables;
}

} // namespace SPLINTER

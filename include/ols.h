/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_OLS_H
#define SPLINTER_OLS_H

#include <definitions.h>
#include <linearfunction.h>
#include <datatable.h>
#include <linearsolvers.h>
#include <utilities.h>

namespace SPLINTER
{

/**
  * Ordinary least-square (OLS)
  */
template <class Vec, class Mat>
DenseVector computeCoefficients(const LinearFunction<Vec, Mat> &func, const DataTable &sample)
{
    // Left hand side
    DenseMatrix X = computeDesignMatrix(func, sample);

    // Right-hand side
    auto y = vectorToDenseVector(sample.getVectorY());

    // Coefficients
    DenseVector c;

    // Solve for coefficients
    DenseQR<> s;
    if (!s.solve(X, y, c))
        throw Exception("computeCoefficients: Failed to solve for coefficients.");

    return c;
};

// TODO: implement OLS with regularization term (lambda/numSample)*coefficients^T*coefficients,
// where lambda >= 0 and default is 0.1/N?
//DenseVector computeCoefficientsRegularized(const LinearFunction &func, const DataTable &samples, double lambda = 0.1);

/**
 * Computes design matrix by evaluating basis functions
 * at each sample point
 */
template<class Vec, class Mat>
Mat computeDesignMatrix(const LinearFunction<Vec, Mat> &func, const DataTable &sample);

template<>
DenseMatrix computeDesignMatrix<DenseVector, DenseMatrix>(const LinearFunction<DenseVector, DenseMatrix> &func, const DataTable &sample)
{
    DenseMatrix X = DenseMatrix::Zero(sample.getNumSamples(), func.getNumCoefficients());

    unsigned int i = 0;
    for (auto it = sample.cbegin(); it != sample.cend(); ++it, ++i)
    {
        // Evaluate basis functions at x
        DenseVector x = vectorToDenseVector(it->getX());
        DenseVector Xi(func.evalBasis(x));

        if (Xi.rows() != func.getNumCoefficients())
            throw Exception("computeDesignMatrix: Xi.rows() != numCoefficients.");

        // Add row to design matrix X
        X.block(i,0,1,func.getNumCoefficients()) = Xi.transpose();
    }

    return X;
};

template<>
SparseMatrix computeDesignMatrix<SparseVector, SparseMatrix>(const LinearFunction<SparseVector, SparseMatrix> &func, const DataTable &sample)
{
    unsigned int numVariables = sample.getNumVariables();
    unsigned int numSamples = sample.getNumSamples();
    unsigned int numCoefficients = func.getNumCoefficients(); // Must equal number of basis functions

    // TODO: Reserve nnz per row (degree+1)
    //int nnzPrCol = bspline.basis.supportedPrInterval();

    SparseMatrix A(numSamples, numCoefficients);
    //A.reserve(DenseVector::Constant(numSamples, nnzPrCol)); // TODO: should reserve nnz per row!

    int i = 0;
    for (auto it = sample.cbegin(); it != sample.cend(); ++it, ++i)
    {
        DenseVector xi(numVariables);
        std::vector<double> xv = it->getX();
        for (unsigned int j = 0; j < numVariables; ++j)
        {
            xi(j) = xv.at(j);
        }

        SparseVector basisValues = func.evalBasis(xi);

        for (SparseVector::InnerIterator it2(basisValues); it2; ++it2)
        {
            A.insert(i,it2.index()) = it2.value();
        }
    }

    A.makeCompressed();

    return A;
};

} // namespace SPLINTER

#endif // SPLINTER_OLS_H

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "ols.h"
#include <linearsolvers.h>

namespace SPLINTER
{

DenseMatrix computeCoefficients(const LinearFunction &func, const DataTable &sample)
{
    // Left hand side
    DenseMatrix X = computeDesignMatrix(func, sample);

    // Right-hand side
    auto yvec = sample.getVectorY();
    DenseVector y(yvec.size());
    for (unsigned int i = 0; i < yvec.size(); ++i)
        y(i) = yvec.at(i);

    // Coefficients
    DenseMatrix c;

    // Solve for coefficients
    DenseQR s;
    if (!s.solve(X, y, c))
        throw Exception("computeCoefficients: Failed to solve for coefficients.");

    // TODO: consider returning a vector instead of a matrix!
    return c;
}

DenseMatrix computeDesignMatrix(const LinearFunction &func, const DataTable &sample)
{
    DenseMatrix X = DenseMatrix::Zero(sample.getNumSamples(), func.getNumCoefficients());

    unsigned int i = 0;
    for (auto it = sample.cbegin(); it != sample.cend(); ++it, ++i)
    {
        // Get point x
        auto xvec = it->getX();
        DenseVector x(xvec.size());
        for (unsigned int j = 0; j < xvec.size(); ++j)
            x(j) = xvec.at(j);

        // Evaluate basis functions at x
        DenseVector Xi(func.evalBasisFunctions(x));

        if (Xi.rows() != func.getNumCoefficients())
        {
            throw Exception("computeDesignMatrix: Xi.rows() != numCoefficients.");
        }

        // Add row to design matrix X
        X.block(i,0,1,func.getNumCoefficients()) = Xi.transpose();
    }

    return X;
}

SparseMatrix computeDesignMatrixSparse(const LinearFunction &func, const DataTable &sample)
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

        SparseVector basisValues = func.evalBasisFunctions(xi);

        for (SparseVector::InnerIterator it2(basisValues); it2; ++it2)
        {
            A.insert(i,it2.index()) = it2.value();
        }
    }

    A.makeCompressed();

    return A;
}

} // namespace SPLINTER

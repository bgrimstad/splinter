/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "ols.h"
//#include <linearsolvers.h>
//#include <utilities.h>

namespace SPLINTER
{

// Full specialization of computeDesignMatrix for dense bases
template<>
DenseMatrix computeDesignMatrix<DenseVector, DenseMatrix>(const LinearFunction<DenseVector, DenseMatrix> &func, const DataTable &data)
{
    DenseMatrix X = DenseMatrix::Zero(data.getNumSamples(), func.getNumCoefficients());

    unsigned int i = 0;
    for (auto it = data.cbegin(); it != data.cend(); ++it, ++i)
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

// Full specialization of computeDesignMatrix for sparse bases
//template<>
//SparseMatrix computeDesignMatrix<SparseVector, SparseMatrix>(const LinearFunction<SparseVector, SparseMatrix> &func, const DataTable &sample)
//{
//    unsigned int numVariables = sample.getNumVariables();
//    unsigned int numSamples = sample.getNumSamples();
//    unsigned int numCoefficients = func.getNumCoefficients(); // Must equal number of basis functions
//
//    // TODO: Reserve nnz per row (degree+1)
//    //int nnzPrCol = bspline.basis.supportedPrInterval();
//
//    SparseMatrix A(numSamples, numCoefficients);
//    //A.reserve(DenseVector::Constant(numSamples, nnzPrCol)); // TODO: should reserve nnz per row!
//
//    int i = 0;
//    for (auto it = sample.cbegin(); it != sample.cend(); ++it, ++i)
//    {
//        DenseVector xi(numVariables);
//        std::vector<double> xv = it->getX();
//        for (unsigned int j = 0; j < numVariables; ++j)
//        {
//            xi(j) = xv.at(j);
//        }
//
//        SparseVector basisValues = func.evalBasis(xi);
//
//        for (SparseVector::InnerIterator it2(basisValues); it2; ++it2)
//        {
//            A.insert(i,it2.index()) = it2.value();
//        }
//    }
//
//    A.makeCompressed();
//
//    return A;
//};

} // namespace SPLINTER

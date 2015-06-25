/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "ordinaryleastsquares.h"
#include "linearsolvers.h"
#include "unsupported/Eigen/KroneckerProduct"

namespace SPLINTER
{

OrdinaryLeastSquares::OrdinaryLeastSquares(const DataTable &samples, unsigned int degree)
    : OrdinaryLeastSquares(samples, std::vector<unsigned int>(samples.getNumVariables(), degree))
{
}

OrdinaryLeastSquares::OrdinaryLeastSquares(const DataTable &samples, std::vector<unsigned int> degrees)
    : degrees(degrees),
      numVariables(samples.getNumVariables()),
      numCoefficients(0)
{
    if (degrees.size() != numVariables)
        throw Exception("OrdinaryLeastSquares::OrdinaryLeastSquares: Inconsistent input data!");

    // Check that a minimum number of samples is provided
    numCoefficients = 1;
    for (auto deg : degrees)
        numCoefficients *= (deg+1);

    if (numCoefficients > samples.getNumSamples())
        throw Exception("OrdinaryLeastSquares::OrdinaryLeastSquares: Insufficient number of samples!");

    // Compute coefficients
    computeCoefficients(samples);
}

double OrdinaryLeastSquares::eval(DenseVector x) const
{
    DenseMatrix res = coefficients*x;
    return res(0,0);
}

void OrdinaryLeastSquares::computeCoefficients(const DataTable &samples)
{
    // Left hand side
    DenseMatrix X = computeDesignMatrix(samples);
    DenseMatrix Xt = X.transpose();
    DenseMatrix XtX = Xt*X;

    // Right-hand side
    auto yvec = samples.getVectorY();
    DenseVector y(yvec.size());
    for (unsigned int i = 0; i < yvec.size(); ++i)
        y(i) = yvec.at(i);
    DenseMatrix Xty = Xt*y;

    // Solve for coefficients
    DenseQR s;
    if (!s.solve(XtX, Xty, coefficients))
        throw Exception("OrdinaryLeastSquares::computeCoefficients: Failed to solve for coefficients.");

    // Transpose coefficients matrix
    coefficients.transposeInPlace();
}

DenseMatrix OrdinaryLeastSquares::computeDesignMatrix(const DataTable &samples) const
{
    DenseMatrix X = DenseMatrix::Zero(samples.getNumSamples(), numCoefficients);

    unsigned int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        // Get point x
        auto xvec = it->getX();
        DenseVector x(xvec.size());
        for (unsigned int j = 0; j < xvec.size(); ++j)
            x(j) = xvec.at(j);

        // Evaluate monomials at x
        DenseVector Xi = evalMonomials(x);
        assert(Xi.cols() == numCoefficients);

        // Add row to design matrix X
        X.block(i,0,1,numCoefficients) = Xi.transpose();
    }

    return X;
}

DenseVector OrdinaryLeastSquares::evalMonomials(DenseVector x) const
{
    std::vector<DenseVector> powers;

    for (unsigned int i = 0; i < numVariables; ++i)
    {
        unsigned int deg = degrees.at(i);
        DenseVector powi = DenseVector::Zero(deg+1);
        for (unsigned int j = 0; j <= deg; ++j)
            powi(j) = std::pow(x(i), j);
        powers.push_back(powi);
    }

    // Kronecker product of monovariable power basis
    DenseVector monomials = DenseVector::Ones(1);

    for (unsigned int i = 0; i < numVariables; ++i)
    {
        DenseVector temp = monomials;
        monomials = kroneckerProduct(temp, powers.at(i));
    }

    assert(monomials.cols() == numCoefficients);

    return monomials;
}

} // namespace SPLINTER

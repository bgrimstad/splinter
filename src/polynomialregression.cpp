/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <serializer.h>
#include "polynomialregression.h"
#include "linearsolvers.h"
#include "mykroneckerproduct.h"

namespace SPLINTER
{

PolynomialRegression::PolynomialRegression()
{
}

PolynomialRegression::PolynomialRegression(const char *fileName)
    : PolynomialRegression(std::string(fileName))
{
}

PolynomialRegression::PolynomialRegression(const std::string fileName)
{
    load(fileName);
}

PolynomialRegression::PolynomialRegression(const DataTable &samples, unsigned int degree)
    : PolynomialRegression(samples, std::vector<unsigned int>(samples.getNumVariables(), degree))
{
}

PolynomialRegression::PolynomialRegression(const DataTable &samples, std::vector<unsigned int> degrees)
    : degrees(degrees),
      numVariables(samples.getNumVariables()),
      numCoefficients(0)
{
    if (degrees.size() != numVariables)
        throw Exception("PolynomialRegression::PolynomialRegression: Inconsistent input data!");

    // Check that a minimum number of samples is provided
    numCoefficients = 1;
    for (auto deg : degrees)
        numCoefficients *= (deg+1);

    if (numCoefficients > samples.getNumSamples())
        throw Exception("PolynomialRegression::PolynomialRegression: Insufficient number of samples!");

    // Compute coefficients
    computeCoefficients(samples);
}

double PolynomialRegression::eval(DenseVector x) const
{
    DenseMatrix monomials = evalMonomials(x);
    DenseMatrix res = coefficients*monomials;
    return res(0,0);
}

DenseMatrix PolynomialRegression::evalJacobian(DenseVector x) const
{
    DenseMatrix jac(1, numVariables);
    jac.fill(0.0);

    for (unsigned int i = 0; i < numVariables; ++i)
    {
        DenseVector diffMonomials = evalDifferentiatedMonomials(x, i);
        DenseVector jaci = coefficients*diffMonomials;
        jac(i) = jaci(0);
    }
    return jac;
}

void PolynomialRegression::computeCoefficients(const DataTable &samples)
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
        throw Exception("PolynomialRegression::computeCoefficients: Failed to solve for coefficients.");

    // Transpose coefficients matrix
    coefficients.transposeInPlace();
}

/*
 * Computes Vandermonde matrix
 * TODO: centering and scaling can improve the numerical properties
 */
DenseMatrix PolynomialRegression::computeDesignMatrix(const DataTable &samples) const
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

        if (Xi.rows() != numCoefficients)
        {
            throw Exception("PolynomialRegression::computeDesignMatrix: Xi.rows() != numCoefficients.");
        }

        // Add row to design matrix X
        X.block(i,0,1,numCoefficients) = Xi.transpose();
    }

    return X;
}

DenseVector PolynomialRegression::evalMonomials(DenseVector x) const
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
    DenseVector monomials = kroneckerProductVectors(powers);

    if (monomials.rows() != numCoefficients)
    {
        throw Exception("PolynomialRegression::evalMonomials: monomials.rows() != numCoefficients.");
    }

    return monomials;
}

DenseVector PolynomialRegression::evalDifferentiatedMonomials(DenseVector x, unsigned int var) const
{
    if (var < 0 || var >= numVariables)
        throw Exception("PolynomialRegression::evalDifferentiatedMonomials: invalid variable.");

    std::vector<DenseVector> powers;

    for (unsigned int i = 0; i < numVariables; ++i)
    {
        unsigned int deg = degrees.at(i);
        DenseVector powi = DenseVector::Zero(deg+1);

        if (var == i)
        {
            // Differentiate wrt. x(i)
            for (unsigned int j = 1; j <= deg; ++j)
                powi(j) = j*std::pow(x(i), j-1);
        }
        else
        {
            for (unsigned int j = 0; j <= deg; ++j)
                powi(j) = std::pow(x(i), j);
        }

        powers.push_back(powi);
    }

    // Kronecker product of monovariable power basis
    DenseVector monomials = kroneckerProductVectors(powers);

    if (monomials.rows() != numCoefficients)
        throw Exception("PolynomialRegression::evalMonomials: monomials.rows() != numCoefficients.");

    return monomials;
}

void PolynomialRegression::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void PolynomialRegression::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string PolynomialRegression::getDescription() const
{
    std::string description("PolynomialRegression of degree");

    // See if all degrees are the same.
    bool equal = true;
    for(size_t i = 1; i < degrees.size(); ++i) {
        equal = equal && (degrees.at(i) == degrees.at(i-1));
    }

    if(equal) {
        description.append(" ");
        description.append(std::to_string(degrees.at(0)));

    } else {
        description.append("s (");
        for(size_t i = 0; i < degrees.size(); ++i) {
            description.append(std::to_string(degrees.at(i)));
            if(i + 1 < degrees.size()) {
                description.append(", ");
            }
        }
        description.append(")");
    }

    return description;
}

} // namespace SPLINTER

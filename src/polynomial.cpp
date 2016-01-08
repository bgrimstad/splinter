/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "polynomial.h"
#include <serializer.h>
#include "mykroneckerproduct.h"
#include <ols.h>

namespace SPLINTER
{

Polynomial::Polynomial(unsigned int numVariables, unsigned int degree)
    : Polynomial(std::vector<unsigned int>(numVariables, degree))
{
}

// Initialize coefficients to zero
Polynomial::Polynomial(std::vector<unsigned int> degrees)
    : LinearFunction<DenseVector, DenseMatrix>(degrees.size(), DenseVector::Zero(computeNumBasisFunctions(degrees))),
      degrees(degrees)
{
}

Polynomial::Polynomial(std::vector<unsigned int> degrees, DenseVector coefficients)
    : LinearFunction<DenseVector, DenseMatrix>(degrees.size(), coefficients),
      degrees(degrees)
{
}

Polynomial::Polynomial(const DataTable &data, unsigned int degree)
        : Polynomial(data, std::vector<unsigned int>(data.getNumVariables(), degree))
{
}

Polynomial::Polynomial(const DataTable &data, std::vector<unsigned int> degrees)
        : Polynomial(degrees)
{
    if (getNumCoefficients() > data.getNumSamples())
        throw Exception("Polynomial::Polynomial: Insufficient number of samples!");

    setCoefficients(computeCoefficients(*this, data));
}

Polynomial::Polynomial(const char *fileName)
        : Polynomial(std::string(fileName))
{
}

Polynomial::Polynomial(const std::string &fileName)
        : LinearFunction(1, DenseVector::Zero(1))
{
    load(fileName);
}

unsigned int Polynomial::computeNumBasisFunctions(std::vector<unsigned int> degrees) const
{
    unsigned int numMonomials = 1;
    for (auto deg : degrees)
        numMonomials *= (deg+1);

    return numMonomials;
}

/**
 * Evaluate monomials
 */
DenseVector Polynomial::evalBasis(DenseVector x) const
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

    if (monomials.rows() != getNumCoefficients())
        throw Exception("Polynomial::evalMonomials: monomials.rows() != numCoefficients.");

    return monomials;
}

DenseMatrix Polynomial::evalBasisJacobian(DenseVector x) const
{
    DenseMatrix jac = DenseMatrix::Zero(getNumCoefficients(), numVariables);

    for (unsigned int var = 0; var < numVariables; ++var)
        jac.block(0,var,getNumCoefficients(),1) = evalDifferentiatedMonomials(x, var);

    return jac;
}

DenseVector Polynomial::evalDifferentiatedMonomials(DenseVector x, unsigned int var) const
{
    if (var >= numVariables)
        throw Exception("Polynomial::evalDifferentiatedMonomials: invalid variable.");

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

    if (monomials.rows() != getNumCoefficients())
        throw Exception("Polynomial::evalMonomials: monomials.rows() != numCoefficients.");

    return monomials;
}

void Polynomial::save(const std::string &fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void Polynomial::load(const std::string &fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

std::string Polynomial::getDescription() const
{
    std::string description("Polynomial of degree");

    // See if all degrees are the same.
    bool equal = true;
    for (size_t i = 1; i < degrees.size(); ++i)
    {
        equal = equal && (degrees.at(i) == degrees.at(i-1));
    }

    if (equal)
    {
        description.append(" ");
        description.append(std::to_string(degrees.at(0)));
    }
    else
    {
        description.append("s (");
        for (size_t i = 0; i < degrees.size(); ++i) {
            description.append(std::to_string(degrees.at(i)));
            if (i + 1 < degrees.size())
            {
                description.append(", ");
            }
        }
        description.append(")");
    }

    return description;
}

} // namespace SPLINTER

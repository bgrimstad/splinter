/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "polynomial2.h"
#include <serializer.h>
#include "mykroneckerproduct.h"
#include <ols.h>

namespace SPLINTER
{

// Initialize coefficients to zero
Polynomial2::Polynomial2(DenseMatrix powers)
    : Polynomial2(powers, DenseVector::Zero(powers.rows()))
{
}

Polynomial2::Polynomial2(DenseMatrix powers, DenseVector coefficients)
    : LinearFunction<DenseVector, DenseMatrix>(powers.cols(), coefficients),
      powers(powers)
{
}

Polynomial2::Polynomial2(const DataTable &data, DenseMatrix powers)
    : Polynomial2(powers)
{
    if (getNumCoefficients() > data.getNumSamples())
        throw Exception("Polynomial::Polynomial: Insufficient number of samples!");

    setCoefficients(computeCoefficients(*this, data));
}

Polynomial2::Polynomial2(const char *fileName)
        : Polynomial2(std::string(fileName))
{
}

Polynomial2::Polynomial2(const std::string &fileName)
        : LinearFunction(1, DenseVector::Zero(1))
{
    load(fileName);
}

/**
 * Evaluate monomials
 */
DenseVector Polynomial2::evalBasis(DenseVector x) const
{
    DenseVector basis = DenseVector::Zero(powers.rows());

    checkInput(x);

    for (unsigned int i = 0; i < powers.rows(); ++i)
    {
        double monomial = 1;
        for (unsigned int j = 0; j < powers.cols(); ++j)
        {
            // Round power to integer
            unsigned int pow = (unsigned int)powers(i, j);
            monomial *= std::pow(x(j), pow);
        }
        basis(i) = monomial;
    }

    return basis;
}

DenseMatrix Polynomial2::evalBasisJacobian(DenseVector x) const
{
    DenseMatrix jac = DenseMatrix::Zero(getNumCoefficients(), numVariables);

    for (unsigned int var = 0; var < numVariables; ++var)
        jac.block(0,var,getNumCoefficients(),1) = evalDifferentiatedMonomials(x, var);

    return jac;
}

DenseVector Polynomial2::evalDifferentiatedMonomials(DenseVector x, unsigned int var) const
{
    if (var >= numVariables)
        throw Exception("Polynomial::evalDifferentiatedMonomials: invalid variable.");

    DenseVector dBasis = DenseVector::Zero(powers.rows());

    checkInput(x);

    for (unsigned int i = 0; i < powers.rows(); ++i)
    {
        double monomial = 1;
        for (unsigned int j = 0; j < powers.cols(); ++j)
        {
            unsigned int pow = (unsigned int)powers(i, j);

            if (var == j)
            {
                if (pow == 0)
                    monomial = 0;
                else
                    monomial *= pow*std::pow(x(j), pow-1);
            }
            else
            {
                monomial *= std::pow(x(j), pow);
            }
        }
        dBasis(i) = monomial;
    }

    return dBasis;
}

unsigned int Polynomial2::getDegree() const
{
    unsigned int maxDeg = 0;

    for (unsigned int i = 0; i < powers.rows(); ++i)
    {
        unsigned int deg = 0;
        for (unsigned int j = 0; j < powers.cols(); ++j)
        {
            deg += powers(i, j);
        }
        if (deg > maxDeg)
            maxDeg = deg;
    }
    return maxDeg;
}

void Polynomial2::save(const std::string &fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void Polynomial2::load(const std::string &fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

std::string Polynomial2::getDescription() const
{
    std::string description("Polynomial = ");

    for (unsigned int i = 0; i < powers.rows(); ++i)
    {
        description.append(" + (");
        description.append(std::to_string(coefficients(i)));
        description.append(")");
        for (unsigned int j = 0; j < powers.cols(); ++j)
        {
            unsigned int pow = powers(i, j);
            if (pow > 0)
            {
                description.append("*x(");
                description.append(std::to_string(j));
                description.append(")^");
                description.append(std::to_string(pow));
            }

        }
    }

    return description;
}

} // namespace SPLINTER

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIAL_H
#define SPLINTER_POLYNOMIAL_H

#include "linearfunction.h"

namespace SPLINTER
{

class SPLINTER_API Polynomial : public LinearFunction
{
public:
    Polynomial(const char *fileName);
    Polynomial(const std::string fileName);
    Polynomial(unsigned int numVariables, unsigned int degree);
    Polynomial(std::vector<unsigned int> degrees);
    Polynomial(std::vector<unsigned int> degrees, DenseVector coefficients);

    SparseVector evalBasisFunctions(DenseVector x) const override;

    SparseMatrix evalBasisFunctionsJacobian(DenseVector x) const override;

    void save(const std::string fileName) const override;

    const std::string getDescription() const;

private:
    Polynomial();

    std::vector<unsigned int> degrees;

    DenseVector evalDifferentiatedMonomials(DenseVector x, unsigned int var) const;

    unsigned int computeNumBasisFunctions(std::vector<unsigned int> degrees) const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const Polynomial &lhs, const Polynomial &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIAL_H

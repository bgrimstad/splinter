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

#include <linearfunction.h>
#include "datatable.h"

namespace SPLINTER
{

class SPLINTER_API Polynomial : public LinearFunction<DenseVector, DenseMatrix>
{
public:
    class Builder;

    Polynomial(std::vector<unsigned int> degrees, DenseVector coefficients);

    /*
     * Construct Polynomial from file
     */
    Polynomial(const char *fileName);
    Polynomial(const std::string &fileName);

    DenseVector evalBasis(DenseVector x) const override;

    DenseMatrix evalBasisJacobian(DenseVector x) const override;

    void save(const std::string &fileName) const override;

    std::string getDescription() const override;

private:
    Polynomial();

    std::vector<unsigned int> degrees;

    DenseVector evalDifferentiatedMonomials(DenseVector x, unsigned int var) const;

    void load(const std::string &fileName) override;

    friend class Serializer;
    friend bool operator==(const Polynomial &lhs, const Polynomial &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIAL_H

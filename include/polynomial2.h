/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIAL2_H
#define SPLINTER_POLYNOMIAL2_H

#include <linearfunction.h>
#include <datatable.h>

namespace SPLINTER
{

/*
 * Class representing polynomials.
 * Powers are given by a matrix P_ij with m rows and n columns,
 * where m equals the number of terms (monomials), and n equals the number of variables.
 * Coefficients are given by a vector of n real numbers.
 * The polynomial can be written as:
 * poly(x) = sum_i (c_i * prod_j x_j^P_ij)
 */
class SPLINTER_API Polynomial2 : public LinearFunction<DenseVector, DenseMatrix>
{
public:
    Polynomial2(DenseMatrix powers);
    Polynomial2(DenseMatrix powers, DenseVector coefficients);
    Polynomial2(const DataTable &data, DenseMatrix powers);
    Polynomial2(const char *fileName);
    Polynomial2(const std::string &fileName);

    DenseVector evalBasis(DenseVector x) const override;

    DenseMatrix evalBasisJacobian(DenseVector x) const override;

    unsigned int getDegree() const;

    void save(const std::string &fileName) const override;

    std::string getDescription() const override;

private:
    Polynomial2();

    DenseMatrix powers;

    DenseVector evalDifferentiatedMonomials(DenseVector x, unsigned int var) const;

    void load(const std::string &fileName) override;

    friend class Serializer;
    friend bool operator==(const Polynomial2 &lhs, const Polynomial2 &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIAL2_H

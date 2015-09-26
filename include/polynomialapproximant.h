/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIALAPPROXIMANT_H
#define SPLINTER_POLYNOMIALAPPROXIMANT_H

#include "approximant.h"
#include "polynomial.h"
#include "datatable.h"

namespace SPLINTER
{

class SPLINTER_API PolynomialApproximant : public Approximant
{
public:
    PolynomialApproximant(const char *fileName);
    PolynomialApproximant(const std::string fileName);
    PolynomialApproximant(const Sample &sample, unsigned int degree);
    PolynomialApproximant(const Sample &sample, std::vector<unsigned int> degrees);

    virtual PolynomialApproximant * clone() const { return new PolynomialApproximant(*this); }

    // Evaluation
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override { DenseMatrix h(numVariables, numVariables); h.fill(0.0); return h; } // TODO: Implement

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

private:
    PolynomialApproximant();

    Polynomial poly;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const PolynomialApproximant &lhs, const PolynomialApproximant &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIALAPPROXIMANT_H

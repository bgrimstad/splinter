/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_POLYNOMIALREGRESSION_H
#define SPLINTER_POLYNOMIALREGRESSION_H

#include "approximant.h"
#include "datatable.h"

namespace SPLINTER
{

class API PolynomialRegression : public Approximant
{
public:
    PolynomialRegression(const char *fileName);
    PolynomialRegression(const std::string fileName);
    PolynomialRegression(const DataTable &samples, unsigned int degree);
    PolynomialRegression(const DataTable &samples, std::vector<unsigned int> degrees);

    virtual PolynomialRegression * clone() const { return new PolynomialRegression(*this); }

    // Evaluation
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override { return DenseMatrix(); } // TODO: Implement
    DenseMatrix evalHessian(DenseVector x) const override { return DenseMatrix(); } // TODO: Implement

    // Getters
    unsigned int getNumVariables() const override { return numVariables; }

    DenseMatrix getCoefficients() const { return coefficients; }

    void save(const std::string fileName) const override;

private:
    PolynomialRegression();

    unsigned int numVariables;
    unsigned int numCoefficients;
    std::vector<unsigned int> degrees;
    DenseMatrix coefficients;

    void computeCoefficients(const DataTable &samples);
    DenseMatrix computeDesignMatrix(const DataTable &samples) const;

    DenseVector evalMonomials(DenseVector x) const;

    void load(const std::string fileName) override;

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_POLYNOMIALREGRESSION_H
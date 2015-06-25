/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_ORDINARYLEASTSQUARES_H
#define SPLINTER_ORDINARYLEASTSQUARES_H

#include "approximant.h"
#include "datatable.h"

namespace SPLINTER
{

class OrdinaryLeastSquares : public Approximant
{
public:
    OrdinaryLeastSquares(const DataTable &samples, unsigned int degree);
    OrdinaryLeastSquares(const DataTable &samples, std::vector<unsigned int> degrees);

    virtual OrdinaryLeastSquares* clone() const { return new OrdinaryLeastSquares(*this); }

    // Evaluation
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override {}
    DenseMatrix evalHessian(DenseVector x) const override {}

    // Getters
    unsigned int getNumVariables() const override { return numVariables; }

    DenseMatrix getCoefficients() const { return coefficients; }

    // Save and load
    void save(const std::string fileName) const override {}
    void load(const std::string fileName) override {}

private:
    unsigned int numVariables;
    unsigned int numCoefficients;
    std::vector<unsigned int> degrees;
    DenseMatrix coefficients;

    void computeCoefficients(const DataTable &samples);
    DenseMatrix computeDesignMatrix(const DataTable &samples) const;

    DenseVector evalMonomials(DenseVector x) const;
};

} // namespace SPLINTER

#endif // SPLINTER_ORDINARYLEASTSQUARES_H

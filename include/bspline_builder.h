/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_BUILDER_H
#define SPLINTER_BSPLINE_BUILDER_H

#include <data_table.h>
#include <knot_utils.h>
#include <bspline.h>

namespace SPLINTER
{

// B-spline smoothing
enum class BSpline::Smoothing {
    NONE,       // No smoothing
    IDENTITY,   // Regularization term alpha*c'*I*c is added to OLS objective
    PSPLINE     // Smoothing term alpha*Delta(c,2) is added to OLS objective
};

// B-spline builder class
class SPLINTER_API BSpline::Builder
{
public:
    Builder(unsigned int dim_x, unsigned int dim_y);

    // Set build options

    Builder& degree(unsigned int degree)
    {
        _degrees = std::vector<unsigned int>(_dim_x, degree);
        return *this;
    }

    Builder& degree(const std::vector<unsigned int> &degrees)
    {
        if (degrees.size() != _dim_x)
            throw Exception("BSpline::Builder::degree: Expected degree vector of length"
                            + std::to_string(_dim_x) + ".");
        _degrees = degrees;
        return *this;
    }

    Builder& numBasisFunctions(unsigned int numBasisFunctions)
    {
        _numBasisFunctions = std::vector<unsigned int>(_dim_x, numBasisFunctions);
        return *this;
    }

    Builder& numBasisFunctions(const std::vector<unsigned int> &numBasisFunctions)
    {
        if (numBasisFunctions.size() != _dim_x)
            throw Exception("BSpline::Builder::numBasisFunctions: Expected numBasisFunctions vector of length "
                            + std::to_string(_dim_x) + ".");
        _numBasisFunctions = numBasisFunctions;
        return *this;
    }

    Builder& knotSpacing(KnotSpacing knotSpacing)
    {
        _knotSpacing = knotSpacing;
        return *this;
    }

    // Fit B-spline to data
//    BSpline fit(const DataTable &data, Smoothing smoothing = Smoothing::NONE, double alpha = .1) const;
    BSpline fit(const DataTable &data,
                Smoothing smoothing = Smoothing::NONE,
                double alpha = .1,
                std::vector<double> weights = std::vector<double>()) const;

private:
    Builder();

    // Control point computations
    DenseMatrix computeControlPoints(const BSpline &bspline,
                                     const DataTable &data,
                                     Smoothing smoothing,
                                     double alpha,
                                     std::vector<double> weights) const;

    SparseMatrix computeBasisFunctionMatrix(const BSpline &bspline, const DataTable &data) const;

    DenseMatrix stackSamplePointValues(const DataTable &data) const;

    // P-spline control point calculation
    SparseMatrix getSecondOrderFiniteDifferenceMatrix(const BSpline &bspline) const;

    // Compute weights matrix from weight vector
    SparseMatrix computeWeightMatrix(const std::vector<double> weights) const;

    // Member variables
    unsigned int _dim_x;
    unsigned int _dim_y;
    std::vector<unsigned int> _degrees;
    std::vector<unsigned int> _numBasisFunctions;
    KnotSpacing _knotSpacing;
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BUILDER_H

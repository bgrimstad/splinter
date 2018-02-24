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
#include <knot_builders.h>
#include <bspline.h>

namespace SPLINTER
{

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

    Builder& num_basis_functions(unsigned int num_basis_functions)
    {
        _numBasisFunctions = std::vector<unsigned int>(_dim_x, num_basis_functions);
        return *this;
    }

    Builder& num_basis_functions(const std::vector<unsigned int> &num_basis_functions)
    {
        if (num_basis_functions.size() != _dim_x)
            throw Exception("BSpline::Builder::num_basis_functions: Expected num_basis_functions vector of length "
                            + std::to_string(_dim_x) + ".");
        _numBasisFunctions = num_basis_functions;
        return *this;
    }

    Builder& knot_spacing(KnotSpacing knot_spacing)
    {
        _knotSpacing = knot_spacing;
        return *this;
    }

    // Fit B-spline to data
    BSpline fit(const DataTable &data,
                Smoothing smoothing = Smoothing::NONE,
                double alpha = .1,
                std::vector<double> weights = std::vector<double>()) const;

private:
    Builder();

    // Member variables
    unsigned int _dim_x;
    unsigned int _dim_y;
    std::vector<unsigned int> _degrees;
    std::vector<unsigned int> _numBasisFunctions;
    KnotSpacing _knotSpacing;
};


/**
 * Convenience functions for B-spline fitting
 */
BSpline bspline_interpolator(const DataTable &data, unsigned int degree = 3);
BSpline bspline_smoother(const DataTable &data, unsigned int degree = 3, double alpha = 0.1);
BSpline pspline_smoother(const DataTable &data, unsigned int degree = 3, double alpha = 0.1);
BSpline cubic_pspline_smoother(const DataTable &data, double alpha = 0.1);


} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BUILDER_H

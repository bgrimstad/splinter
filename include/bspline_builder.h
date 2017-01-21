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

#include "data_table.h"
#include "bspline.h"

namespace SPLINTER
{

// B-spline smoothing
enum class BSpline::Smoothing
{
    NONE,       // No smoothing
    IDENTITY,   // Regularization term alpha*c'*I*c is added to OLS objective
    PSPLINE     // Smoothing term alpha*Delta(c,2) is added to OLS objective
};

/*
 * B-spline knot spacing
 */
enum class BSpline::KnotSpacing
{
    /*
     * Clamped and with knots that mimic the spacing of sample points using a moving average filter.
     * p+1 multiplicity of end knots.
     * Example:
     *      Given sample points [0, 1, 2, 3, 4, 5]
     *      Then, for degree 3 the knot vector becomes: [0, 0, 0, 0, 2, 3, 5, 5, 5, 5]
     *      and for degree 1 the knot vector becomes: [0, 0, 1, 2, 3, 4, 5, 5]
     */
    AS_SAMPLED,

    /*
     * Clamped knot vector with equidistant internal knots. p+1 multiplicity of end knots.
     * Example:
     *      Given samples on the interval [0, 5]
     *      For degree 3: [0, 0, 0, 0, 1, 2, 3, 4, 5, 5, 5, 5]
     *      For degree 1: [0, 0, 1, 2, 3, 4, 5, 5]
     */
    EQUIDISTANT,

    /*
     * Experimental knot spacing for testing purposes only.
     * Currently, it gives a non-clamped knot vector of equidistant knots.
     * NOTE: may eventually become EQUIDISTANT_NOT_CLAMPED
     * Example:
     *      For any degree: [0, 1, 2, 3, 4, 5]
     *
     */
    EXPERIMENTAL
};

// B-spline builder class
class SPLINTER_API BSpline::Builder
{
public:
    Builder(unsigned int dim_x, unsigned int dim_y);

    // Set build options

    Builder& degree(unsigned int degree)
    {
        _degrees = getBSplineDegrees(_dim_x, degree);
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
    BSpline fit(const DataTable &data, Smoothing smoothing = Smoothing::NONE, double alpha = .1) const;

private:
    Builder();

    std::vector<unsigned int> getBSplineDegrees(unsigned int numVariables, unsigned int degree)
    {
        // TODO: Remove this test
        if (degree > 5)
            throw Exception("BSpline::Builder: Only degrees in range [0, 5] are supported.");
        return std::vector<unsigned int>(numVariables, degree);
    }

    // Control point computations
    DenseMatrix computeControlPoints(const BSpline &bspline,
                                     const DataTable &data,
                                     Smoothing smoothing,
                                     double alpha) const;

    SparseMatrix computeBasisFunctionMatrix(const BSpline &bspline, const DataTable &data) const;

    DenseMatrix stackSamplePointValues(const DataTable &data) const;

    // P-spline control point calculation
    SparseMatrix getSecondOrderFiniteDifferenceMatrix(const BSpline &bspline) const;

    // Computing knots
    std::vector<std::vector<double>> computeKnotVectors(const DataTable &data) const;

    std::vector<double> computeKnotVector(const std::vector<double> &values, unsigned int degree,
                                          unsigned int numBasisFunctions) const;

    // Member variables
    unsigned int _dim_x;
    unsigned int _dim_y;
    std::vector<unsigned int> _degrees;
    std::vector<unsigned int> _numBasisFunctions;
    KnotSpacing _knotSpacing;
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BUILDER_H

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
    Builder(const DataTable &data);

    Builder& alpha(double alpha)
    {
        if (alpha < 0)
            throw Exception("BSpline::Builder::alpha: alpha must be non-negative.");

        _alpha = alpha;
        return *this;
    }

    // Set build options

    Builder& degree(unsigned int degree)
    {
        _degrees = getBSplineDegrees(_data.getDimX(), degree);
        return *this;
    }

    Builder& degree(const std::vector<unsigned int> &degrees)
    {
        if (degrees.size() != _data.getDimX())
            throw Exception("BSpline::Builder: Inconsistent length on degree vector.");
        _degrees = degrees;
        return *this;
    }

    Builder& numBasisFunctions(unsigned int numBasisFunctions)
    {
        _numBasisFunctions = std::vector<unsigned int>(_data.getDimX(), numBasisFunctions);
        return *this;
    }

    Builder& numBasisFunctions(const std::vector<unsigned int> &numBasisFunctions)
    {
        if (numBasisFunctions.size() != _data.getDimX())
            throw Exception("BSpline::Builder: Inconsistent length on numBasisFunctions vector.");
        _numBasisFunctions = numBasisFunctions;
        return *this;
    }

    Builder& knotSpacing(KnotSpacing knotSpacing)
    {
        _knotSpacing = knotSpacing;
        return *this;
    }

    Builder& smoothing(Smoothing smoothing)
    {
        _smoothing = smoothing;
        return *this;
    }

    // Build B-spline
    BSpline build() const;

private:
    Builder();

    std::vector<unsigned int> getBSplineDegrees(unsigned int numVariables, unsigned int degree)
    {
        if (degree > 5)
            throw Exception("BSpline::Builder: Only degrees in range [0, 5] are supported.");
        return std::vector<unsigned int>(numVariables, degree);
    }

    // Control point computations
    DenseMatrix computeControlPoints(const BSpline &bspline) const;
    SparseMatrix computeBasisFunctionMatrix(const BSpline &bspline) const;
    DenseMatrix stackSamplePointValues() const;
    // P-spline control point calculation
    SparseMatrix getSecondOrderFiniteDifferenceMatrix(const BSpline &bspline) const;

    // Computing knots
    std::vector<std::vector<double>> computeKnotVectors() const;
    std::vector<double> computeKnotVector(const std::vector<double> &values, unsigned int degree, unsigned int numBasisFunctions) const;

    // Member variables
    DataTable _data;
    std::vector<unsigned int> _degrees;
    std::vector<unsigned int> _numBasisFunctions;
    KnotSpacing _knotSpacing;
    Smoothing _smoothing;
    double _alpha;
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BUILDER_H

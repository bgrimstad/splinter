/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEBUILDER_H
#define SPLINTER_BSPLINEBUILDER_H

#include "datatable.h"
#include "bspline.h"

namespace SPLINTER
{
    // B-spline degrees
    enum class BSpline::Degree
    {
        LINEAR,     // Linear basis functions in each variable
        QUADRATIC,  // Quadratic basis functions in each variable
        CUBIC,      // Cubic basis functions in each variable
        QUARTIC     // Quartic basis functions in each variable
    };

    // B-spline smoothing
    enum class BSpline::Smoothing
    {
        NONE,           // No smoothing
        REGULARIZATION, // Regularization term lambda*c^2 is added to OLS objective
        PSPLINE         // Smoothing term lambda*Delta(c,2) is added to OLS objective
    };

    // B-spline knot spacing
    enum class BSpline::KnotSpacing
    {
        SAMPLE,         // Knot spacing mimicking sample spacing (moving average)
        EQUIDISTANT,    // Equidistant knots
        EXPERIMENTAL    // Experimental knot spacing (needs more testing)
    };

    // B-spline builder class
    class SPLINTER_API BSpline::Builder
    {
    public:
        Builder(const DataTable &data);

        // Set build options
        Builder& degree(Degree degree)
        {
            _degrees = getBSplineDegrees(_data.getNumVariables(), degree);
            return *this;
        }

        Builder& degree(unsigned int numVariables, Degree degree)
        {
            _degrees = getBSplineDegrees(numVariables, degree);
            return *this;
        }

        Builder& numKnots(unsigned int numVariables, unsigned int numKnots)
        {
            _numKnots = std::vector<unsigned int>(numVariables, numKnots);
            return *this;
        }

        Builder& numKnots(std::vector<unsigned int> numKnots)
        {
            _numKnots = numKnots;
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

        Builder& lambda(double lambda)
        {
            if (lambda < 0)
                throw Exception("BSpline::Builder::lambda: Lambda must be non-negative.");

            _lambda = lambda;
            return *this;
        }

        // Build B-spline
        BSpline build() const;

    private:
        Builder();

        std::vector<unsigned int> getBSplineDegrees(unsigned int numVariables, Degree degree)
        {
            switch (degree)
            {
                case Degree::LINEAR:
                    return std::vector<unsigned int>(numVariables, 1);
                case Degree::QUADRATIC:
                    return std::vector<unsigned int>(numVariables, 2);
                case Degree::CUBIC:
                    return std::vector<unsigned int>(numVariables, 3);
                case Degree::QUARTIC:
                    return std::vector<unsigned int>(numVariables, 4);
                default:
                    return std::vector<unsigned int>(numVariables, 3);
            }
        }

        // Control point computations
        virtual DenseMatrix computeCoefficients(const DataTable &samples, const BSpline &bspline) const;
        virtual DenseMatrix computeBSplineCoefficients(const DataTable &samples, const BSpline &bspline) const;
        SparseMatrix computeBasisFunctionMatrix(const DataTable &samples, const BSpline &bspline) const;
        DenseMatrix controlPointEquationRHS(const DataTable &samples) const;

        // P-spline control point calculation
        DenseMatrix computePSplineCoefficients(const DataTable &samples, const BSpline &bspline) const;
        SparseMatrix getSecondOrderFiniteDifferenceMatrix(const BSpline &bspline) const;

        // Computing knots
        std::vector<std::vector<double> > computeKnotVectors(const DataTable &data, std::vector<unsigned int> degrees) const;
        std::vector<double> computeKnotVector(const std::vector<double> &values, unsigned int degree) const;
        std::vector<double> knotVectorMovingAverage(const std::vector<double> &values, unsigned int degree) const;
        std::vector<double> knotVectorEquidistant(const std::vector<double> &values, unsigned int degree) const;
        std::vector<double> knotVectorBuckets(const std::vector<double> &values, unsigned int degree, unsigned int maxSegments = 10) const;

        // Auxiliary
        std::vector<double> extractUniqueSorted(const std::vector<double> &values) const;

        // Member variables
        DataTable _data;
        std::vector<unsigned int> _degrees;
        std::vector<unsigned int> _numKnots;
        KnotSpacing _knotSpacing;
        Smoothing _smoothing;
        double _lambda;
    };

} // namespace SPLINTER

#endif //SPLINTER_BSPLINEBUILDER_H

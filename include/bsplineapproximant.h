/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEAPPROXIMANT_H
#define SPLINTER_BSPLINEAPPROXIMANT_H

#include "datatable.h"
#include "bspline.h"

namespace SPLINTER
{

// Enum for different B-spline types
enum class BSplineType
{
    LINEAR,     // Linear basis functions in each variable
    QUADRATIC,  // Quadratic basis functions in each variable
    CUBIC,      // Cubic basis functions in each variable
    QUARTIC     // Quartic basis functions in each variable
};

inline std::vector<unsigned int> getBSplineDegrees(unsigned int numVariables, BSplineType type)
{
    switch (type) {
    case BSplineType::LINEAR:
            return std::vector<unsigned int>(numVariables, 1);
        break;
    case BSplineType::QUADRATIC:
            return std::vector<unsigned int>(numVariables, 2);
        break;
    case BSplineType::CUBIC:
            return std::vector<unsigned int>(numVariables, 3);
        break;
    case BSplineType::QUARTIC:
            return std::vector<unsigned int>(numVariables, 4);
        break;
    default:
        // Default is CUBIC
        return std::vector<unsigned int>(numVariables, 3);
        break;
    }
}

class SPLINTER_API BSplineApproximant : public Approximant
{
public:
    BSplineApproximant(const char *fileName);
    BSplineApproximant(const std::string fileName);
    BSplineApproximant(const DataTable &samples, std::vector<unsigned int> basisDegrees);
    BSplineApproximant(const DataTable &samples, BSplineType type);

    virtual BSplineApproximant * clone() const { return new BSplineApproximant(*this); }

    // Build B-spline
    BSpline buildBSpline(const DataTable &samples, std::vector<unsigned int> basisDegrees) const;

    // Evaluation
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override;

    BSpline getBSpline() const
    {
        return bspline;
    }

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

protected:
    BSplineApproximant();
    BSplineApproximant(unsigned int numVariables);

    BSpline bspline;

    // Control point computations
    virtual DenseMatrix computeCoefficients(const DataTable &samples) const;
    SparseMatrix computeBasisFunctionMatrix(const DataTable &samples) const;
    DenseMatrix controlPointEquationRHS(const DataTable &samples) const;

    // Computing knots
    std::vector<std::vector<double> > computeKnotVectorsFromSamples(const DataTable &samples, std::vector<unsigned int> degrees) const;
    virtual std::vector<double> computeKnotVector(const std::vector<double> &values, unsigned int degree) const;
    std::vector<double> knotVectorMovingAverage(const std::vector<double> &values, unsigned int degree) const;
    std::vector<double> knotVectorEquidistant(const std::vector<double> &values, unsigned int degree) const;
    std::vector<double> knotVectorBuckets(const std::vector<double> &values, unsigned int degree, unsigned int maxSegments = 10) const;

    // Auxiliary
    std::vector<double> extractUniqueSorted(const std::vector<double> &values) const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const BSplineApproximant &lhs, const BSplineApproximant &rhs);
};


} // namespace SPLINTER

#endif // SPLINTER_BSPLINEAPPROXIMANT_H

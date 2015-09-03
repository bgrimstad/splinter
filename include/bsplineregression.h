/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINEREGRESSION_H
#define SPLINTER_BSPLINEREGRESSION_H

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

/**
 * Class that implements the multivariate tensor product B-spline
 */
class SPLINTER_API BSplineRegression : public Approximant
{
public:

    BSplineRegression(unsigned int numVaribles);

    /**
     * Construct B-spline that interpolates the samples in DataTable
     */
    BSplineRegression(const DataTable &samples, std::vector<unsigned int> basisDegrees);
    BSplineRegression(const DataTable &samples, BSplineType type);

    /**
     * Construct B-spline from file
     */
    BSplineRegression(const char *fileName);
    BSplineRegression(const std::string fileName);

    virtual BSplineRegression* clone() const { return new BSplineRegression(*this); }

    // Evaluation
    double eval(DenseVector x) const override;
    DenseMatrix evalJacobian(DenseVector x) const override;
    DenseMatrix evalHessian(DenseVector x) const override;

    void save(const std::string fileName) const override;

    // Getters
    BSpline getBSpline() const
    {
        return bspline;
    }

    const std::string getDescription() const override;

protected:
    BSplineRegression();

    BSpline bspline;

    // Control point computations
    virtual DenseMatrix computeControlPoints(const DataTable &samples);
    SparseMatrix computeBasisFunctionMatrix(const DataTable &samples) const;
    void controlPointEquationRHS(const DataTable &samples, DenseMatrix &Bx, DenseMatrix &By) const;

    // Computing knots
    std::vector<std::vector<double> > computeKnotVectorsFromSamples(const DataTable &samples, std::vector<unsigned int> degrees) const;
    std::vector<double> knotVectorMovingAverage(std::vector<double> &vec, unsigned int degree) const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const BSplineRegression &lhs, const BSplineRegression &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINEREGRESSION_H

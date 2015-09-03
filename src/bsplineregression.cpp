/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bsplineregression.h"
#include "mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "linearsolvers.h"
#include <serializer.h>
#include <iostream>

namespace SPLINTER
{

BSplineRegression::BSplineRegression()
    : Approximant(1),
      bspline(BSpline(1))
{
}

BSplineRegression::BSplineRegression(const char *fileName)
    : BSplineRegression(std::string(fileName))
{
}

BSplineRegression::BSplineRegression(const std::string fileName)
    : Approximant(1),
      bspline(BSpline(1))
{
    load(fileName);
}

BSplineRegression::BSplineRegression(const DataTable &samples, std::vector<unsigned int> basisDegrees)
    : Approximant(samples.getNumVariables()),
      bspline(BSpline(computeKnotVectorsFromSamples(samples, basisDegrees), basisDegrees))
{
    // Check data
    if (!samples.isGridComplete())
        throw Exception("BSpline::BSpline: Cannot create B-spline from irregular (incomplete) grid.");

    bspline.setCoefficients(computeControlPoints(samples));
}

BSplineRegression::BSplineRegression(const DataTable &samples, BSplineType type = BSplineType::CUBIC)
    : Approximant(samples.getNumVariables()),
      bspline(BSpline(samples.getNumVariables()))
{
    // Check data
    if (!samples.isGridComplete())
        throw Exception("BSpline::BSpline: Cannot create B-spline from irregular (incomplete) grid.");

    // Default is CUBIC
    std::vector<unsigned int> basisDegrees(samples.getNumVariables(), 3);

    // Set multivariate basis
    if (type == BSplineType::LINEAR)
        basisDegrees = std::vector<unsigned int>(samples.getNumVariables(), 1);
    else if (type == BSplineType::QUADRATIC)
        basisDegrees = std::vector<unsigned int>(samples.getNumVariables(), 2);
    else if (type == BSplineType::CUBIC)
        basisDegrees = std::vector<unsigned int>(samples.getNumVariables(), 3);
    else if (type == BSplineType::QUARTIC)
        basisDegrees = std::vector<unsigned int>(samples.getNumVariables(), 4);

    // Compute knot vectors from samples
    auto knotVectors = computeKnotVectorsFromSamples(samples, basisDegrees);

    // Create B-spline
    bspline = BSpline(knotVectors, basisDegrees);

    // Compute control points
    bspline.setCoefficients(computeControlPoints(samples));
}

double BSplineRegression::eval(DenseVector x) const
{
    return bspline.eval(x);
}

DenseMatrix BSplineRegression::evalJacobian(DenseVector x) const
{
    return bspline.evalJacobian(x);
}

DenseMatrix BSplineRegression::evalHessian(DenseVector x) const
{
    return bspline.evalHessian(x);
}

DenseMatrix BSplineRegression::computeControlPoints(const DataTable &samples)
{
    /* Setup and solve equations Ac = b,
     * A = basis functions at sample x-values,
     * b = sample y-values when calculating control coefficients,
     * b = sample x-values when calculating knot averages
     * c = control coefficients or knot averages.
     */
    SparseMatrix A = computeBasisFunctionMatrix(samples);

    DenseMatrix Bx, By;
    controlPointEquationRHS(samples, Bx, By);

    DenseMatrix Cx, Cy;

    int numEquations = A.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if (!solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using sparse solver." << std::endl;
#endif // NDEBUG

        SparseLU s;
        bool successfulSolve = (s.solve(A,Bx,Cx) && s.solve(A,By,Cy));

        solveAsDense = !successfulSolve;
    }

    if (solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using dense solver." << std::endl;
#endif // NDEBUG

        DenseMatrix Ad = A.toDense();
        DenseQR s;
        bool successfulSolve = (s.solve(Ad,Bx,Cx) && s.solve(Ad,By,Cy));
        if (!successfulSolve)
        {
            throw Exception("BSpline::computeControlPoints: Failed to solve for B-spline coefficients.");
        }
    }

    return Cy.transpose();
    //knotaverages = Cx.transpose();
}

SparseMatrix BSplineRegression::computeBasisFunctionMatrix(const DataTable &samples) const
{
    unsigned int numVariables = samples.getNumVariables();
    unsigned int numSamples = samples.getNumSamples();

    // TODO: Reserve nnz per row (degree+1)
    //int nnzPrCol = bspline.basis.supportedPrInterval();

    SparseMatrix A(numSamples, bspline.getNumBasisFunctionsTotal());
    //A.reserve(DenseVector::Constant(numSamples, nnzPrCol)); // TODO: should reserve nnz per row!

    int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        DenseVector xi(numVariables);
        std::vector<double> xv = it->getX();
        for (unsigned int j = 0; j < numVariables; ++j)
        {
            xi(j) = xv.at(j);
        }

        SparseVector basisValues = bspline.evalBasisFunctions(xi);

        for (SparseVector::InnerIterator it2(basisValues); it2; ++it2)
        {
            A.insert(i,it2.index()) = it2.value();
        }
    }

    A.makeCompressed();

    return A;
}

void BSplineRegression::controlPointEquationRHS(const DataTable &samples, DenseMatrix &Bx, DenseMatrix &By) const
{
    unsigned int numVariables = samples.getNumVariables();
    unsigned int numSamples = samples.getNumSamples();

    Bx.resize(numSamples, numVariables);
    By.resize(numSamples, 1);

    int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        std::vector<double> x = it->getX();

        for (unsigned int j = 0; j < x.size(); ++j)
        {
            Bx(i,j) = x.at(j);
        }

        By(i,0) = it->getY();
    }
}

std::vector<std::vector<double> > BSplineRegression::computeKnotVectorsFromSamples(const DataTable &samples, std::vector<unsigned int> degrees) const
{
    if (samples.getNumVariables() != degrees.size())
        throw Exception("BSpline::computeKnotVectorsFromSamples: Inconsistent sizes on input vectors.");

    std::vector<std::vector<double> > grid = samples.getTableX();

    std::vector<std::vector<double> > knotVectors;

    for (unsigned int i = 0; i < samples.getNumVariables(); ++i)
    {
        // Using moving average filter to compute knot vectors
        auto knotVec = knotVectorMovingAverage(grid.at(i), degrees.at(i));

        knotVectors.push_back(knotVec);
    }

    return knotVectors;
}

/*
 * Automatic construction of (p+1)-regular knot vector
 * using moving average.
 *
 * Requirement:
 * Knot vector should be of size n+p+1.
 * End knots are should be repeated p+1 times.
 *
 * Computed sizes:
 * n+2*(p) = n + p + 1 + (p - 1)
 * k = (p - 1) values must be removed from sample vector.
 * w = k + 3 window size in moving average
 *
 * Algorithm:
 * 1) compute n - k values using moving average with window size w
 * 2) repeat first and last value p + 1 times
 *
 * The resulting knot vector has n - k + 2*p = n + p + 1 knots.
 *
 * NOTE:
 * For _equidistant_ samples, the resulting knot vector is identicaly to
 * the free end conditions knot vector used in cubic interpolation.
 * That is, samples (a,b,c,d,e,f) produces the knot vector (a,a,a,a,c,d,f,f,f,f) for p = 3.
 * For p = 1, (a,b,c,d,e,f) becomes (a,a,b,c,d,e,f,f).
 *
 */
std::vector<double> BSplineRegression::knotVectorMovingAverage(std::vector<double> &vec, unsigned int degree) const
{
    // Sort and remove duplicates
    std::vector<double> uniqueX(vec);
    std::sort(uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy(uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize(distance(uniqueX.begin(),it));

    // Compute sizes
    unsigned int n = uniqueX.size();
    unsigned int k = degree-1; // knots to remove
    unsigned int w = k + 3; // Window size

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "BSpline::knotVectorMovingAverage: Only " << n
          << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
          << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    std::vector<double> knots(n-k-2, 0);

    // Compute (n-k-2) interior knots using moving average
    for (unsigned int i = 0; i < n-k-2; ++i)
    {
        double ma = 0;
        for (unsigned int j = 0; j < w; ++j)
            ma += uniqueX.at(i+j);

        knots.at(i) = ma/w;
    }

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.begin(), uniqueX.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.end(), uniqueX.back());

    // Number of knots in a (p+1)-regular knot vector
    assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

void BSplineRegression::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void BSplineRegression::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string BSplineRegression::getDescription() const
{
    return bspline.getDescription();
}

} // namespace SPLINTER

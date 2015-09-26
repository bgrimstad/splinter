/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bsplineapproximant.h"
#include "mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"
#include <linearsolvers.h>
#include <serializer.h>
#include <iostream>
#include <utilities.h>

namespace SPLINTER
{

BSplineApproximant::BSplineApproximant(unsigned int numVariables)
    : Approximant(numVariables),
      bspline(numVariables)
{
}

BSplineApproximant::BSplineApproximant(const Sample &samples, std::vector<unsigned int> basisDegrees)
    : Approximant(samples.getNumVariables()),
      bspline(buildBSpline(samples, basisDegrees))
{
    // Compute coefficients
    auto coefficients = computeCoefficients(samples);
    bspline.setCoefficients(coefficients);
}

// TODO: this may build the B-spline twice!
BSplineApproximant::BSplineApproximant(const Sample &samples, BSplineType type = BSplineType::CUBIC)
    : BSplineApproximant(samples, getBSplineDegrees(samples.getNumVariables(), type))
{
}

/*
 * Construct from saved data
 */
BSplineApproximant::BSplineApproximant(const char *fileName)
    : BSplineApproximant(std::string(fileName))
{
}

BSplineApproximant::BSplineApproximant(const std::string fileName)
    : Approximant(1),
      bspline(1)
{
    load(fileName);
}

/*
 * Build B-spline
 */
BSpline BSplineApproximant::buildBSpline(const Sample &samples, std::vector<unsigned int> basisDegrees) const
{
    // Check data
    if (!samples.isGridComplete())
        throw Exception("BSplineApproximant::buildBSpline: Cannot create B-spline from irregular (incomplete) grid.");

    auto knotVectors = computeKnotVectorsFromSamples(samples, basisDegrees);

    return BSpline(knotVectors, basisDegrees);
}

double BSplineApproximant::eval(DenseVector x) const
{
    return bspline.eval(x);
}

DenseMatrix BSplineApproximant::evalJacobian(DenseVector x) const
{
    return bspline.evalJacobian(x);
}

DenseMatrix BSplineApproximant::evalHessian(DenseVector x) const
{
    return bspline.evalHessian(x);
}

DenseMatrix BSplineApproximant::computeCoefficients(const Sample &samples) const
{
    /* Setup and solve equations Ac = b,
     * A = basis functions at sample x-values,
     * b = sample y-values when calculating control coefficients,
     * b = sample x-values when calculating knot averages
     * c = control coefficients or knot averages.
     */
    SparseMatrix A2 = computeBasisFunctionMatrix(samples);
    DenseMatrix b2 = controlPointEquationRHS(samples);

    SparseMatrix At = A2.transpose();
    SparseMatrix A = At*A2; // Multiply with transpose to obtain a symmetric matrix
    DenseMatrix b = At*b2;

    DenseMatrix w;

    int numEquations = A.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    // TODO: compute only coefficients (knot averages not needed)
    if (!solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using sparse solver." << std::endl;
#endif // NDEBUG

        SparseLU s;
        //bool successfulSolve = (s.solve(A,Bx,Cx) && s.solve(A,By,Cy));

        solveAsDense = !s.solve(A,b,w);
    }

    if (solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using dense solver." << std::endl;
#endif // NDEBUG

        DenseMatrix Ad = A.toDense();
        DenseQR s;
        //bool successfulSolve = (s.solve(Ad,Bx,Cx) && s.solve(Ad,By,Cy));
        if (!s.solve(Ad,b,w))
        {
            throw Exception("BSpline::computeControlPoints: Failed to solve for B-spline coefficients.");
        }
    }

    return w.transpose();
}

SparseMatrix BSplineApproximant::computeBasisFunctionMatrix(const Sample &samples) const
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

DenseMatrix BSplineApproximant::controlPointEquationRHS(const Sample &samples) const
{
    DenseMatrix B = DenseMatrix::Zero(samples.getNumSamples(), 1);

    int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
        B(i,0) = it->getY();

    return B;
}

std::vector<std::vector<double> > BSplineApproximant::computeKnotVectorsFromSamples(const Sample &samples, std::vector<unsigned int> degrees) const
{
    if (samples.getNumVariables() != degrees.size())
        throw Exception("BSpline::computeKnotVectorsFromSamples: Inconsistent sizes on input vectors.");

    std::vector<std::vector<double> > grid = samples.getTableX();

    std::vector<std::vector<double> > knotVectors;

    for (unsigned int i = 0; i < samples.getNumVariables(); ++i)
    {
        // Using moving average filter to compute knot vectors
        auto knotVec = computeKnotVector(grid.at(i), degrees.at(i));

        knotVectors.push_back(knotVec);
    }

    return knotVectors;
}

std::vector<double> BSplineApproximant::computeKnotVector(const std::vector<double> &values, unsigned int degree) const
{
    return knotVectorMovingAverage(values, degree);
    //return knotVectorBuckets(values, degree);
    //return knotVectorEquidistant(values, degree);
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
 * TODO: does not work well when number of knots is << number of samples! For such cases
 * almost all knots will lie close to the left samples. Try a bucket approach, where the
 * samples are added to buckets and the knots computed as the average of these.
 */
std::vector<double> BSplineApproximant::knotVectorMovingAverage(const std::vector<double> &values, unsigned int degree) const
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // Compute sizes
    unsigned int n = unique.size();
    unsigned int k = degree-1; // knots to remove
    unsigned int w = k + 3; // Window size

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knotVectorMovingAverage: Only " << n
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
            ma += unique.at(i+j);

        knots.at(i) = ma/w;
    }

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

std::vector<double> BSplineApproximant::knotVectorEquidistant(const std::vector<double> &values, unsigned int degree) const
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // Compute sizes
    unsigned int n = unique.size();
    unsigned int k = degree-1; // knots to remove

    // The minimum number of samples from which a free knot vector can be created
    if (n < degree+1)
    {
        std::ostringstream e;
        e << "knotVectorMovingAverage: Only " << n
          << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
          << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    // Compute (n-k-2) equidistant interior knots
    unsigned int numIntKnots = std::max(n-k-2, (unsigned int)0);
    numIntKnots = std::min((unsigned int)10, numIntKnots);
    std::vector<double> knots = linspace(unique.front(), unique.back(), numIntKnots);

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (unsigned int i = 0; i < degree; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree; ++i)
        knots.insert(knots.end(), unique.back());

    // Number of knots in a (p+1)-regular knot vector
    //assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

std::vector<double> BSplineApproximant::knotVectorBuckets(const std::vector<double> &values, unsigned int degree, unsigned int maxSegments) const
{
    // Sort and remove duplicates
    std::vector<double> unique = extractUniqueSorted(values);

    // The minimum number of samples from which a free knot vector can be created
    if (unique.size() < degree+1)
    {
        std::ostringstream e;
        e << "BSplineApproximant::knotVectorBuckets: Only " << unique.size()
          << " unique sample points are given. A minimum of degree+1 = " << degree+1
          << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    // Num internal knots (0 <= ni <= unique.size() - degree - 1)
    unsigned int ni = unique.size() - degree - 1;

    // Num segments
    unsigned int ns = ni + degree + 1;

    // Limit number of segments
    if (ns > maxSegments && maxSegments >= degree + 1)
    {
        ns = maxSegments;
        ni = ns - degree - 1;
    }

    // Num knots
    unsigned int nk = ns + degree + 1;

    // Check numbers
    if (ni < 0 || ni > unique.size() - degree - 1)
        throw Exception("BSplineApproximant::knotVectorBuckets: Invalid number of internal knots!");

    // Compute window sizes
    unsigned int w = 0;
    if (ni > 0)
        w = std::floor(unique.size()/ni);

    // Residual
    unsigned int res = unique.size() - w*ni;

    // Create array with window sizes
    std::vector<unsigned int> windows(ni, w);

    // Add residual
    for (unsigned int i = 0; i < res; ++i)
        windows.at(i) += 1;

    // Compute internal knots
    std::vector<double> knots(ni, 0);

    // Compute (n-k-2) interior knots using moving average
    unsigned int index = 0;
    for (unsigned int i = 0; i < ni; ++i)
    {
        for (unsigned int j = 0; j < windows.at(i); ++j)
        {
            knots.at(i) += unique.at(index+j);
        }
        knots.at(i) /= windows.at(i);
        index += windows.at(i);
    }

    // Repeat first knot p + 1 times (for interpolation of start point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.begin(), unique.front());

    // Repeat last knot p + 1 times (for interpolation of end point)
    for (unsigned int i = 0; i < degree + 1; ++i)
        knots.insert(knots.end(), unique.back());

    return knots;
}

std::vector<double> BSplineApproximant::extractUniqueSorted(const std::vector<double> &values) const
{
    // Sort and remove duplicates
    std::vector<double> unique(values);
    std::sort(unique.begin(), unique.end());
    std::vector<double>::iterator it = unique_copy(unique.begin(), unique.end(), unique.begin());
    unique.resize(distance(unique.begin(),it));
    return unique;
}

void BSplineApproximant::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void BSplineApproximant::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string BSplineApproximant::getDescription() const
{
    std::string description("BSplineApproximant");

    description.append(bspline.getDescription());

    return description;
}

} // namespace SPLINTER

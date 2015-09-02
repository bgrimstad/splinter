/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline.h"
#include "bsplinebasis.h"
#include "mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "linearsolvers.h"
#include <serializer.h>
#include <iostream>

namespace SPLINTER
{

BSpline::BSpline()
    : Approximant(1)
{}

/*
 * Constructors for multivariate B-splines with explicit data
 */
BSpline::BSpline(std::vector<double> coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees)
    : Approximant(knotVectors.size()),
      coefficients(DenseMatrix::Zero(1, coefficients.size()))
{
    for (unsigned int i = 0; i < coefficients.size(); ++i)
        this->coefficients(0,i) = coefficients[i];

    if (this->coefficients.rows() != 1)
        throw Exception("BSpline::BSpline: coefficient matrix can only have one row!");

    basis = BSplineBasis(knotVectors, basisDegrees);

    computeKnotAverages();

    init();

    checkControlPoints();
}

BSpline::BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees)
    : Approximant(knotVectors.size()),
      coefficients(coefficients)
{
    if (coefficients.rows() != 1)
        throw Exception("BSpline::BSpline: coefficient matrix can only have one row!");

    basis = BSplineBasis(knotVectors, basisDegrees);

    computeKnotAverages();

    init();

    checkControlPoints();
}

/*
 * Constructors for interpolation of samples in DataTable
 */
BSpline::BSpline(const DataTable &samples, unsigned int degree)
    : Approximant(samples.getNumVariables())
{
    // Check data
    if (!samples.isGridComplete())
        throw Exception("BSpline::BSpline: Cannot create B-spline from irregular (incomplete) grid.");

    // Assuming that all basis function are of the same degree
    std::vector<unsigned int> basisDegrees(samples.getNumVariables(), degree);

    // Compute knot vectors from samples
    auto knotVectors = computeKnotVectorsFromSamples(samples, basisDegrees);

    // Set multivariate basis
    basis = BSplineBasis(knotVectors, basisDegrees);

    // Calculate control points
    computeControlPoints(samples);

    init();

    checkControlPoints();
}

BSpline::BSpline(const DataTable &samples, BSplineType type = BSplineType::CUBIC)
    : Approximant(samples.getNumVariables())
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

    // Set multivariate basis
    basis = BSplineBasis(knotVectors, basisDegrees);

    // Calculate control points
    computeControlPoints(samples);

    init();

    checkControlPoints();
}

/*
 * Construct from saved data
 */
BSpline::BSpline(const char *fileName)
    : BSpline(std::string(fileName))
{
}

BSpline::BSpline(const std::string fileName)
    : Approximant(1)
{
    load(fileName);
}

void BSpline::init()
{
    bool initialKnotRefinement = false;
    if (initialKnotRefinement)
    {
        globalKnotRefinement();
        checkControlPoints();
    }
}

double BSpline::eval(DenseVector x) const
{
	if (!pointInDomain(x))
	{
		throw Exception("BSpline::eval: Evaluation at point outside domain.");
	}

	SparseVector tensorvalues = basis.eval(x);
	DenseVector y = coefficients*tensorvalues;
	return y(0);
}

double BSpline::eval(double x) const
{
	DenseVector x_vec(1);
	x_vec(0) = x;
	return this->eval(x_vec);
}

/*
 * Returns the Jacobian evaluated at x.
 * The Jacobian is an 1 x n matrix,
 * where n is the dimension of x.
 */
DenseMatrix BSpline::evalJacobian(DenseVector x) const
{
    if (!pointInDomain(x))
    {
        throw Exception("BSpline::evalJacobian: Evaluation at point outside domain.");
    }

    //SparseMatrix Bi = basis.evalBasisJacobian(x);       // Sparse Jacobian implementation
    //SparseMatrix Bi = basis.evalBasisJacobian2(x);  // Sparse Jacobian implementation
    DenseMatrix Bi = basis.evalBasisJacobianOld(x);  // Old Jacobian implementation

    return coefficients*Bi;
}

/*
 * Returns the Hessian evaluated at x.
 * The Hessian is an n x n matrix,
 * where n is the dimension of x.
 */
DenseMatrix BSpline::evalHessian(DenseVector x) const
{
    if (!pointInDomain(x))
    {        
        throw Exception("BSpline::evalHessian: Evaluation at point outside domain.");
    }

    DenseMatrix H;
    H.setZero(1,1);
    DenseMatrix identity = DenseMatrix::Identity(numVariables,numVariables);
    DenseMatrix caug;
    caug = kroneckerProduct(identity, coefficients);
    DenseMatrix DB = basis.evalBasisHessian(x);
    H = caug*DB;

    // Fill in upper triangular of Hessian
    for (size_t i = 0; i < numVariables; ++i)
    {
        for (size_t j = i+1; j < numVariables; ++j)
        {
            H(i,j) = H(j,i);
        }
    }
    return H;
}

std::vector<unsigned int> BSpline::getNumBasisFunctions() const
{
    std::vector<unsigned int> ret;
    for (unsigned int i = 0; i < numVariables; i++)
        ret.push_back(basis.getNumBasisFunctions(i));
    return ret;
}

std::vector< std::vector<double> > BSpline::getKnotVectors() const
{
    return basis.getKnotVectors();
}

std::vector<unsigned int> BSpline::getBasisDegrees() const
{
    return basis.getBasisDegrees();
}

std::vector<double> BSpline::getDomainUpperBound() const
{
    return basis.getSupportUpperBound();
}

std::vector<double> BSpline::getDomainLowerBound() const
{
    return basis.getSupportLowerBound();
}

DenseMatrix BSpline::getControlPoints() const
{
    int nc = coefficients.cols();
    DenseMatrix controlPoints(numVariables + 1, nc);

    controlPoints.block(0, 0, numVariables, nc) = knotaverages;
    controlPoints.block(numVariables, 0, 1, nc) = coefficients;

    return controlPoints;
}

void BSpline::setControlPoints(DenseMatrix &controlPoints)
{
    assert(controlPoints.rows() == numVariables + 1);
    int nc = controlPoints.cols();

    knotaverages = controlPoints.block(0, 0, numVariables, nc);
    coefficients = controlPoints.block(numVariables, 0, 1, nc);

    checkControlPoints();
}

void BSpline::checkControlPoints() const
{
    if (coefficients.cols() != knotaverages.cols())
        throw Exception("BSpline::checkControlPoints: Inconsistent size of coefficients and knot averages matrices.");
    if (knotaverages.rows() != numVariables)
        throw Exception("BSpline::checkControlPoints: Inconsistent size of knot averages matrix.");
    if (coefficients.rows() != 1)
        throw Exception("BSpline::checkControlPoints: Coefficients matrix does not have one row.");
}

bool BSpline::pointInDomain(DenseVector x) const
{
    return basis.insideSupport(x);
}

void BSpline::reduceDomain(std::vector<double> lb, std::vector<double> ub, bool doRegularizeKnotVectors)
{
    if (lb.size() != numVariables || ub.size() != numVariables)
        throw Exception("BSpline::reduceDomain: Inconsistent vector sizes!");

    std::vector<double> sl = basis.getSupportLowerBound();
    std::vector<double> su = basis.getSupportUpperBound();

    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        // Check if new domain is empty
        if (ub.at(dim) <= lb.at(dim) || lb.at(dim) >= su.at(dim) || ub.at(dim) <= sl.at(dim))
            throw Exception("BSpline::reduceDomain: Cannot reduce B-spline domain to empty set!");

        // Check if new domain is a strict subset
        if (su.at(dim) < ub.at(dim) || sl.at(dim) > lb.at(dim))
            throw Exception("BSpline::reduceDomain: Cannot expand B-spline domain!");

        // Tightest possible
        sl.at(dim) = lb.at(dim);
        su.at(dim) = ub.at(dim);
    }

    if (doRegularizeKnotVectors)
    {
        regularizeKnotVectors(sl, su);
    }

    // Remove knots and control points that are unsupported with the new bounds
    if (!removeUnsupportedBasisFunctions(sl, su))
    {
        throw Exception("BSpline::reduceDomain: Failed to remove unsupported basis functions!");
    }
}

void BSpline::globalKnotRefinement()
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.refineKnots();

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();
}

void BSpline::localKnotRefinement(DenseVector x)
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.refineKnotsLocally(x);

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();
}

void BSpline::decomposeToBezierForm()
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.decomposeToBezierForm();

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();
}

// Computes knot averages: assumes that basis is initialized!
void BSpline::computeKnotAverages()
{
    // Calculate knot averages for each knot vector
    std::vector< DenseVector > knotAverages;
    for (unsigned int i = 0; i < numVariables; i++)
    {
        std::vector<double> knots = basis.getKnotVector(i);
        DenseVector mu; mu.setZero(basis.getNumBasisFunctions(i));

        for (unsigned int j = 0; j < basis.getNumBasisFunctions(i); j++)
        {
            double knotAvg = 0;
            for (unsigned int k = j+1; k <= j+basis.getBasisDegree(i); k++)
            {
                knotAvg += knots.at(k);
            }
            mu(j) = knotAvg/basis.getBasisDegree(i);
        }
        knotAverages.push_back(mu);
    }

    // Calculate vectors of ones (with same length as corresponding knot average vector)
    std::vector<DenseVector> knotOnes;
    for (unsigned int i = 0; i < numVariables; i++)
    {
        DenseVector ones;
        ones.setOnes(knotAverages.at(i).rows());
        knotOnes.push_back(ones);
    }

    // Fill knot average matrix one row at the time
    // NOTE: Must have same pattern as samples in DataTable
    // TODO: Use DataTable to achieve the same pattern
    knotaverages.resize(numVariables, basis.getNumBasisFunctions());

    for (unsigned int i = 0; i < numVariables; i++)
    {
        DenseMatrix mu_ext(1,1); mu_ext(0,0) = 1;
        for (unsigned int j = 0; j < numVariables; j++)
        {
            DenseMatrix temp = mu_ext;
            if (i == j)
                mu_ext = Eigen::kroneckerProduct(temp, knotAverages.at(j));
            else
                mu_ext = Eigen::kroneckerProduct(temp, knotOnes.at(j));
        }
        assert(mu_ext.rows() == basis.getNumBasisFunctions());
        knotaverages.block(i,0,1,basis.getNumBasisFunctions()) = mu_ext.transpose();
    }


    assert(knotaverages.rows() == numVariables && knotaverages.cols() == basis.getNumBasisFunctions());
}

void BSpline::computeControlPoints(const DataTable &samples)
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

    coefficients = Cy.transpose();
    knotaverages = Cx.transpose();
}

SparseMatrix BSpline::computeBasisFunctionMatrix(const DataTable &samples) const
{
    unsigned int numVariables = samples.getNumVariables();
    unsigned int numSamples = samples.getNumSamples();

    int nnzPrCol = basis.supportedPrInterval();

    SparseMatrix A(numSamples, basis.getNumBasisFunctions());
    A.reserve(DenseVector::Constant(numSamples, nnzPrCol)); // TODO: should reserve nnz per row!

    int i = 0;
    for (auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        DenseVector xi(numVariables);
        std::vector<double> xv = it->getX();
        for (unsigned int j = 0; j < numVariables; ++j)
        {
            xi(j) = xv.at(j);
        }

        SparseVector basisValues = basis.eval(xi);

        for (SparseVector::InnerIterator it2(basisValues); it2; ++it2)
        {
            A.insert(i,it2.index()) = it2.value();
        }
    }

    A.makeCompressed();

    return A;
}

void BSpline::controlPointEquationRHS(const DataTable &samples, DenseMatrix &Bx, DenseMatrix &By) const
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

void BSpline::insertKnots(double tau, unsigned int dim, unsigned int multiplicity)
{
    // Insert knots and compute knot insertion matrix
    SparseMatrix A = basis.insertKnots(tau, dim, multiplicity);

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();
}

void BSpline::regularizeKnotVectors(std::vector<double> &lb, std::vector<double> &ub)
{
    // Add and remove controlpoints and knots to make the b-spline p-regular with support [lb, ub]
    if (!(lb.size() == numVariables && ub.size() == numVariables))
        throw Exception("BSpline::regularizeKnotVectors: Inconsistent vector sizes.");

    for (unsigned int dim = 0; dim < numVariables; dim++)
    {
        unsigned int multiplicityTarget = basis.getBasisDegree(dim) + 1;

        // Inserting many knots at the time (to save number of B-spline coefficient calculations)
        // NOTE: This method generates knot insertion matrices with more nonzero elements than
        // the method that inserts one knot at the time. This causes the preallocation of
        // kronecker product matrices to become too small and the speed deteriorates drastically
        // in higher dimensions because reallocation is necessary. This can be prevented by
        // precomputing the number of nonzeros when preallocating memory (see myKroneckerProduct).
        int numKnotsLB = multiplicityTarget - basis.getKnotMultiplicity(dim, lb.at(dim));
        if (numKnotsLB > 0)
        {
            insertKnots(lb.at(dim), dim, numKnotsLB);
        }

        int numKnotsUB = multiplicityTarget - basis.getKnotMultiplicity(dim, ub.at(dim));
        if (numKnotsUB > 0)
        {
            insertKnots(ub.at(dim), dim, numKnotsUB);
        }
    }
}

bool BSpline::removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub)
{
    assert(lb.size() == numVariables);
    assert(ub.size() == numVariables);

    SparseMatrix A = basis.reduceSupport(lb, ub);

    if (coefficients.cols() != A.rows())
        return false;

    // Remove unsupported control points (basis functions)
    coefficients = coefficients*A;
    knotaverages = knotaverages*A;

    return true;
}

void BSpline::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void BSpline::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string BSpline::getDescription() const
{
    std::string description("BSpline of degree");
    auto degrees = getBasisDegrees();
    // See if all degrees are the same.
    bool equal = true;
    for (size_t i = 1; i < degrees.size(); ++i)
    {
        equal = equal && (degrees.at(i) == degrees.at(i-1));
    }

    if(equal)
    {
        description.append(" ");
        description.append(std::to_string(degrees.at(0)));
    }
    else
    {
        description.append("s (");
        for (size_t i = 0; i < degrees.size(); ++i)
        {
            description.append(std::to_string(degrees.at(i)));
            if (i + 1 < degrees.size())
            {
                description.append(", ");
            }
        }
        description.append(")");
    }

    return description;
}

std::vector<std::vector<double> > BSpline::computeKnotVectorsFromSamples(const DataTable &samples, std::vector<unsigned int> degrees) const
{
    if (samples.getNumVariables() != degrees.size())
        throw Exception("BSpline::computeKnotVectorsFromSamples: Inconsistent sizes on input vectors.");

    std::vector< std::vector<double> > grid = samples.getTableX();

    std::vector<std::vector<double> > knotVectors;

    for (unsigned int i = 0; i < numVariables; ++i)
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
std::vector<double> BSpline::knotVectorMovingAverage(std::vector<double> &vec, unsigned int degree) const
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

} // namespace SPLINTER

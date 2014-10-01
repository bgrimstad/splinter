/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "bspline.h"
#include "include/bsplinebasis.h"
#include "include/mykroneckerproduct.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "include/linearsolvers.h"

#include <iostream>

namespace MultivariateSplines
{

// Constructor for explicitly given multivariate B-splines
BSpline::BSpline(DenseMatrix coefficients, std::vector< std::vector<double> > knotVectors, std::vector<unsigned int> basisDegrees)
    : coefficients(coefficients)
{
    numVariables = knotVectors.size();
    assert(coefficients.rows() == 1);

    basis = BSplineBasis(knotVectors, basisDegrees, KnotVectorType::EXPLICIT);

    computeKnotAverages();

    init();

    checkControlPoints();
}

// Constructors for interpolation of samples in DataTable
BSpline::BSpline(const DataTable &samples, BSplineType type = BSplineType::CUBIC_FREE)
{
    // Check data
    if(!samples.isGridComplete())
    {
        throw Exception("BSpline::BSpline: Cannot create B-spline from irregular (incomplete) grid.");
    }

    numVariables = samples.getNumVariables();

    std::vector< std::vector<double> > xdata = samples.getTableX();

    // Set multivariate basis
    if(type == BSplineType::LINEAR)
    {
        std::vector<unsigned int> basisDegrees(samples.getNumVariables(), 1);
        basis = BSplineBasis(xdata, basisDegrees, KnotVectorType::FREE);
    }
    else if(type == BSplineType::QUADRATIC_FREE)
    {
        std::vector<unsigned int> basisDegrees(samples.getNumVariables(), 2);
        basis = BSplineBasis(xdata, basisDegrees, KnotVectorType::FREE);
    }
    else if(type == BSplineType::CUBIC_FREE)
    {
        std::vector<unsigned int> basisDegrees(samples.getNumVariables(), 3);
        basis = BSplineBasis(xdata, basisDegrees, KnotVectorType::FREE);
    }
    else
    {
        // Default is CUBIC_FREE
        std::vector<unsigned int> basisDegrees(samples.getNumVariables(), 3);
        basis = BSplineBasis(xdata, basisDegrees, KnotVectorType::FREE);
    }

    // Calculate control points
    computeControlPoints(samples);

    init();

    checkControlPoints();
}

void BSpline::init()
{
    bool initialKnotRefinement = false;
    if(initialKnotRefinement)
    {
        refineKnotVectors();
        checkControlPoints();
    }
}

double BSpline::eval(DenseVector x) const
{
    if(!pointInDomain(x))
    {
        throw Exception("BSpline::eval: Evaluation at point outside domain.");
    }

    SparseVector tensorvalues = basis.eval(x);
    DenseVector y = coefficients*tensorvalues;
    return y(0);
}

/*
 * Returns the Jacobian evaluated at x.
 * The Jacobian is an 1 x n matrix,
 * where n is the dimension of x.
 */
DenseMatrix BSpline::evalJacobian(DenseVector x) const
{
    if(!pointInDomain(x))
    {
        throw Exception("BSpline::evalJacobian: Evaluation at point outside domain.");
    }

    DenseMatrix Bi = basis.evalBasisJacobianOld(x); // Old Jacobian implementation
    //SparseMatrix Ji = basis.evalBasisJacobian(x);   // New Jacobian implementation

//    // Test difference in Jacobians
//    DenseMatrix Jid(Ji);
//    DenseMatrix dJ = Bi - Jid;
//    DenseVector errorVec = dJ.rowwise().maxCoeff();
//    DenseVector error = errorVec.colwise().maxCoeff();
//    if (std::abs(error(0)) > 1e-10)
//    {
//        std::cout << "NOTABLE DIFFERENCE IN JACOBIANS: " << std::abs(error(0)) << std::endl;
//        exit(1);
//    }

    return coefficients*Bi;
}

/*
 * Returns the Hessian evaluated at x.
 * The Hessian is an n x n matrix,
 * where n is the dimension of x.
 */
DenseMatrix BSpline::evalHessian(DenseVector x) const
{
    if(!pointInDomain(x))
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
    return H;
}

std::vector< std::vector<double> > BSpline::getKnotVectors() const
{
    return basis.getKnotVectors();
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

bool BSpline::checkControlPoints() const
{
    assert(coefficients.cols() == knotaverages.cols());
    assert(knotaverages.rows() == numVariables);
    assert(coefficients.rows() == 1);
    return true;
}

bool BSpline::pointInDomain(DenseVector x) const
{
    return basis.insideSupport(x);
}

bool BSpline::reduceDomain(std::vector<double> lb, std::vector<double> ub, bool doRegularizeKnotVectors, bool doRefineKnotVectors)
{
    if(lb.size() != numVariables || ub.size() != numVariables)
        return false;

    std::vector<double> sl = basis.getSupportLowerBound();
    std::vector<double> su = basis.getSupportUpperBound();

    bool isStrictSubset = false;

    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        if(ub.at(dim) <= lb.at(dim)
            || lb.at(dim) >= su.at(dim)
            || ub.at(dim) <= sl.at(dim))
        {
            throw Exception("BSpline::reduceDomain: Cannot reduce B-spline domain to empty set!");
        }

        if(su.at(dim) > ub.at(dim))
        {
            isStrictSubset = true;
            su.at(dim) = ub.at(dim);
        }

        if(lb.at(dim) > sl.at(dim))
        {
            isStrictSubset = true;
            sl.at(dim) = lb.at(dim);
        }
    }

    if(isStrictSubset)
    {
        if (doRegularizeKnotVectors && !regularizeKnotVectors(sl, su))
        {
            throw Exception("BSpline::reduceDomain: Failed to regularize knot vectors!");
        }

        // Remove knots and control points that are unsupported with the new bounds
        if(!removeUnsupportedBasisFunctions(sl, su))
        {
            throw Exception("BSpline::reduceDomain: Failed to remove unsupported basis functions!");
        }

        // Refine knots
        if(doRefineKnotVectors && !refineKnotVectors())
        {
            throw Exception("BSpline::reduceDomain: Failed to refine knot vectors!");
        }
    }

    return true;
}

// Computes knot averages: assumes that basis is initialized!
void BSpline::computeKnotAverages()
{
    // Calculate knot averages for each knot vector
    std::vector< DenseVector > knotAverages;
    for(unsigned int i = 0; i < numVariables; i++)
    {
        std::vector<double> knots = basis.getKnotVector(i);
        DenseVector mu; mu.setZero(basis.numBasisFunctions(i));

        for (unsigned int j = 0; j < basis.numBasisFunctions(i); j++)
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
    for(unsigned int i = 0; i < numVariables; i++)
    {
        DenseVector ones;
        ones.setOnes(knotAverages.at(i).rows());
        knotOnes.push_back(ones);
    }

    // Fill knot average matrix one row at the time
    // NOTE: Must have same pattern as samples in DataTable
    // TODO: Use DataTable to achieve the same pattern
    knotaverages.resize(numVariables, basis.numBasisFunctions());

    for(unsigned int i = 0; i < numVariables; i++)
    {
        DenseMatrix mu_ext(1,1); mu_ext(0,0) = 1;
        for (unsigned int j = 0; j < numVariables; j++)
        {
            DenseMatrix temp = mu_ext;
            if(i == j)
                mu_ext = Eigen::kroneckerProduct(temp, knotAverages.at(j));
            else
                mu_ext = Eigen::kroneckerProduct(temp, knotOnes.at(j));
        }
        assert(mu_ext.rows() == basis.numBasisFunctions());
        knotaverages.block(i,0,1,basis.numBasisFunctions()) = mu_ext.transpose();
    }


    assert(knotaverages.rows() == numVariables && knotaverages.cols() == basis.numBasisFunctions());
}

void BSpline::computeControlPoints(const DataTable &samples)
{
    /* Setup and solve equations Ac = b,
     * A = basis functions at sample x-values,
     * b = sample y-values when calculating control coefficients,
     * b = sample x-values when calculating knot averages
     * c = control coefficients or knot averages.
     */
    SparseMatrix A;
    computeBasisFunctionMatrix(samples, A);

    DenseMatrix Bx, By;
    controlPointEquationRHS(samples, Bx, By);

    DenseMatrix Cx, Cy;

    int numEquations = A.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if(!solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using sparse solver." << std::endl;
#endif // NDEBUG

        SparseLU s;
        bool successfulSolve = (s.solve(A,Bx,Cx) && s.solve(A,By,Cy));

        solveAsDense = !successfulSolve;
    }

    if(solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using dense solver." << std::endl;
#endif // NDEBUG

        DenseMatrix Ad = A.toDense();
        DenseQR s;
        bool successfulSolve = (s.solve(Ad,Bx,Cx) && s.solve(Ad,By,Cy));
        if(!successfulSolve)
        {
            throw Exception("BSpline::computeControlPoints: Failed to solve for B-spline coefficients.");
        }
    }

    coefficients = Cy.transpose();
    knotaverages = Cx.transpose();
}

void BSpline::computeBasisFunctionMatrix(const DataTable &samples, SparseMatrix &A) const
{
    unsigned int numVariables = samples.getNumVariables();
    unsigned int numSamples = samples.getNumSamples();

    int nnzPrCol = basis.supportedPrInterval();

    A.resize(numSamples, basis.numBasisFunctions());
    A.reserve(DenseVector::Constant(numSamples, nnzPrCol)); // TODO: should reserve nnz per row!

    int i = 0;
    for(auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        DenseVector xi(numVariables);
        std::vector<double> xv = it->getX();
        for(unsigned int j = 0; j < numVariables; ++j)
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
}

void BSpline::controlPointEquationRHS(const DataTable &samples, DenseMatrix &Bx, DenseMatrix &By) const
{
    unsigned int numVariables = samples.getNumVariables();
    unsigned int numSamples = samples.getNumSamples();

    Bx.resize(numSamples, numVariables);
    By.resize(numSamples, 1);

    int i = 0;
    for(auto it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        std::vector<double> x = it->getX();

        for(unsigned int j = 0; j < x.size(); ++j)
        {
            Bx(i,j) = x.at(j);
        }

        By(i,0) = it->getY();
    }
}

bool BSpline::insertKnots(double tau, unsigned int dim, unsigned int multiplicity)
{
    // Test multiplicity at knot
    if(basis.getKnotMultiplicity(dim, tau) + multiplicity > basis.getBasisDegree(dim) + 1)
        return false;

    // Insert knots and compute knot insertion matrix
    SparseMatrix A;
    if(!basis.insertKnots(A, tau, dim, multiplicity))
        return false;

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();

    return true;
}

bool BSpline::refineKnotVectors()
{
    // Compute knot insertion matrix
    SparseMatrix A;
    if(!basis.refineKnots(A))
        return false;

    // Update control points
    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();

    return true;
}

// NOTE: Do this lower in the hierarchy so that all knots can be added at once
bool BSpline::regularizeKnotVectors(std::vector<double> &lb, std::vector<double> &ub)
{
    // Add and remove controlpoints and knots to make the b-spline p-regular with support [lb, ub]
    if(!(lb.size() == numVariables && ub.size() == numVariables))
        return false;

    for(unsigned int dim = 0; dim < numVariables; dim++)
    {
        unsigned int multiplicityTarget = basis.getBasisDegree(dim) + 1;

        // Inserting many knots at the time (to save number of B-spline coefficient calculations)
        // NOTE: This method generates knot insertion matrices with more nonzero elements than
        // the method that inserts one knot at the time. This causes the preallocation of
        // kronecker product matrices to become too small and the speed deteriorates drastically
        // in higher dimensions because reallocation is necessary. This can be prevented by
        // precomputing the number of nonzeros when preallocating memory (see myKroneckerProduct).
        int numKnotsLB = multiplicityTarget - basis.getKnotMultiplicity(dim, lb.at(dim));
        if(numKnotsLB > 0)
        {
            if(!insertKnots(lb.at(dim), dim, numKnotsLB))
                return false;
        }

        int numKnotsUB = multiplicityTarget - basis.getKnotMultiplicity(dim, ub.at(dim));
        if(numKnotsUB > 0)
        {
            if(!insertKnots(ub.at(dim), dim, numKnotsUB))
                return false;
        }

        // Old insertion method: inserts one knot at the time
//        while (basis.getKnotMultiplicity(dim, lb.at(dim)) < multiplicityTarget)
//        {
//            insertKnots(lb.at(dim), dim);
//        }

//        while (basis.getKnotMultiplicity(dim, ub.at(dim)) < multiplicityTarget)
//        {
//            insertKnots(ub.at(dim), dim);
//        }
    }

    return true;
}

bool BSpline::removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub)
{
    assert(lb.size() == numVariables);
    assert(ub.size() == numVariables);

    SparseMatrix A;
    if(!basis.reduceSupport(lb, ub, A))
        return false;

    if(coefficients.cols() != A.rows())
        return false;

    // Remove unsupported control points (basis functions)
    coefficients = coefficients*A;
    knotaverages = knotaverages*A;

    return true;
}

} // namespace MultivariateSplines

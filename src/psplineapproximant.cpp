/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <serializer.h>
#include "psplineapproximant.h"
#include "linearsolvers.h"

namespace SPLINTER
{

PSplineApproximant::PSplineApproximant(const Sample &samples, std::vector<unsigned int> basisDegrees, double lambda)
    : BSplineApproximant(samples, basisDegrees),
      lambda(lambda)
{
    // Check lambda
    if (lambda <= 0)
        throw Exception("buildPSpline: Lambda must be strictly positive.");

    // Compute coefficients
    // TODO: coefficients are computed twice since BSplineApproximant constructor calls BSplineApproximant::computesCoefficients
    auto coefficients = computeCoefficients(samples);
    bspline.setCoefficients(coefficients);
}

PSplineApproximant::PSplineApproximant(const Sample &samples, BSplineType type, double lambda)
    : PSplineApproximant(samples, getBSplineDegrees(samples.getNumVariables(), type), lambda)
{
}

PSplineApproximant::PSplineApproximant(const Sample &samples, double lambda)
    : PSplineApproximant(samples, std::vector<unsigned int>(samples.getNumVariables(), 3), lambda)
{
}

/*
 * Construct from saved data
 */
PSplineApproximant::PSplineApproximant(const char *fileName)
    : PSplineApproximant(std::string(fileName))
{
}

PSplineApproximant::PSplineApproximant(const std::string fileName)
    : BSplineApproximant(1),
      lambda(0.03)
{
    load(fileName);
}

DenseMatrix PSplineApproximant::computeCoefficients(const Sample &samples) const
{
    // Assuming regular grid
    unsigned int numSamples = samples.getNumSamples();

    /* Setup and solve equations Lc = R,
     * L = B'*W*B + l*D'*D
     * R = B'*W*y
     * c = control coefficients or knot averages.
     * B = basis functions at sample x-values,
     * W = weighting matrix for interpolating specific points
     * D = second-order finite difference matrix
     * l = penalizing parameter (increase for more smoothing)
     * y = sample y-values when calculating control coefficients,
     * y = sample x-values when calculating knot averages
     */

    SparseMatrix L, W;

    // Weight matrix
    W.resize(numSamples, numSamples);
    W.setIdentity();

    // Basis function matrix
    SparseMatrix B = computeBasisFunctionMatrix(samples);

    // Second order finite difference matrix
    SparseMatrix D = getSecondOrderFiniteDifferenceMatrix();

    // Left-hand side matrix
    L = B.transpose()*W*B + lambda*D.transpose()*D;

    // Compute right-hand side matrices
    DenseMatrix By = controlPointEquationRHS(samples);
    //Rx = B.transpose()*W*Bx;
    DenseMatrix Ry = B.transpose()*W*By;

    // Matrices to store the resulting coefficients
    DenseMatrix Cy;

    int numEquations = L.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if (!solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using sparse solver." << std::endl;
#endif // NDEBUG

        SparseLU s;
        bool successfulSolve = s.solve(L,Ry,Cy);

        solveAsDense = !successfulSolve;
    }

    if (solveAsDense)
    {
#ifndef NDEBUG
        std::cout << "Computing B-spline control points using dense solver." << std::endl;
#endif // NDEBUG

        DenseMatrix Ld = L.toDense();
        DenseQR s;
        bool successfulSolve = s.solve(Ld, Ry, Cy);

        if (!successfulSolve)
        {
            throw Exception("PSpline::computeControlPoints: Failed to solve for B-spline coefficients.");
        }
    }

    return Cy.transpose();
}

/*
 * Function for generating second order finite-difference matrix, which is used for penalizing the
 * (approximate) second derivative in control point calculation for P-splines.
 */
SparseMatrix PSplineApproximant::getSecondOrderFiniteDifferenceMatrix() const
{
    unsigned int numVariables = bspline.getNumVariables();

    // Number of (total) basis functions - defines the number of columns in D
    unsigned int numCols = bspline.getNumBasisFunctionsTotal();
    std::vector<unsigned int> numBasisFunctions = bspline.getNumBasisFunctions();

    // Number of basis functions (and coefficients) in each variable
    std::vector<unsigned int> dims;
    for (unsigned int i = 0; i < numVariables; i++)
        dims.push_back(numBasisFunctions.at(i));

    std::reverse(dims.begin(), dims.end());

    for (unsigned int i=0; i < numVariables; i++)
    {
        // Need at least three coefficients in each variable
        assert(basis.getNumBasisFunctions(i) >= 3);
    }

    // Number of rows in D and in each block
    int numRows = 0;
    std::vector< int > numBlkRows;
    for (unsigned int i = 0; i < numVariables; i++)
    {
        int prod = 1;
        for (unsigned int j = 0; j < numVariables; j++)
        {
            if (i == j)
                prod *= (dims[j] - 2);
            else
                prod *= dims[j];
        }
        numRows += prod;
        numBlkRows.push_back(prod);
    }

    // Resize and initialize D
    SparseMatrix D(numRows, numCols);
    D.reserve(DenseVector::Constant(numCols,2*numVariables));   // D has no more than two elems per col per dim

    int i = 0;                                          // Row index
    // Loop though each dimension (each dimension has its own block)
    for (unsigned int d = 0; d < numVariables; d++)
    {
        // Calculate left and right products
        int leftProd = 1;
        int rightProd = 1;
        for (unsigned int k = 0; k < d; k++)
        {
            leftProd *= dims[k];
        }
        for (unsigned int k = d+1; k < numVariables; k++)
        {
            rightProd *= dims[k];
        }

        // Loop through subblocks on the block diagonal
        for (int j = 0; j < rightProd; j++)
        {
            // Start column of current subblock
            int blkBaseCol = j*leftProd*dims[d];
            // Block rows [I -2I I] of subblock
            for (unsigned int l = 0; l < (dims[d] - 2); l++)
            {
                // Special case for first dimension
                if (d == 0)
                {
                    int k = j*leftProd*dims[d] + l;
                    D.insert(i,k) = 1;
                    k += leftProd;
                    D.insert(i,k) = -2;
                    k += leftProd;
                    D.insert(i,k) = 1;
                    i++;
                }
                else
                {
                    // Loop for identity matrix
                    for (int n = 0; n < leftProd; n++)
                    {
                        int k = blkBaseCol + l*leftProd + n;
                        D.insert(i,k) = 1;
                        k += leftProd;
                        D.insert(i,k) = -2;
                        k += leftProd;
                        D.insert(i,k) = 1;
                        i++;
                    }
                }
            }
        }
    }

//    // Number of (total) basis functions - defines the number of columns in D
//    int numCols = basis.getNumBasisFunctions();

//    // Number of basis functions (and coefficients) in each dimension
//    std::vector < int > dims = basis.getTensorIndexDimension();
//    std::reverse(dims.begin(), dims.end()); // flip vector

//    // Number of variables
//    int vars = dims.size();

//    for (int i=0; i < vars; i++)
//    {
//        // Need at least three coefficients in each dimension
//        assert(basis.getNumBasisFunctions(i) >= 3);
//        dims.push_back(basis.getNumBasisFunctions(i));
//    }

//    // Number of rows in D and in each block
//    int numRows = 0;
//    std::vector< int > numBlkRows;
//    for (int i = 0; i < vars; i++)
//    {
//        int prod = 1;
//        for (int j = 0; j < vars; j++)
//        {
//            if (i == j)
//                prod *= (dims[j] - 2);
//            else
//                prod *= dims[j];
//        }
//        numRows += prod;
//        numBlkRows.push_back(prod);
//    }

//    // Resize and initialize D
//    DenseMatrix Dd;
//    Dd.resize(numRows, numCols);
//    Dd.setZero();
//    cout << numRows << "/" << numCols << endl;

//    // Accumulate variable for product of preceding dimensions
//    int acc = 1;

//    // Accumulate variable for block insertion point
//    int insertionRow = 0;

//    // Base block
//    DenseMatrix d(1,3);
//    d << 1, -2, 1;

//    // Calculate block matrix associated with each variable and insert into D
//    for (int i = 0; i < vars; i++)
//    {

//        DenseMatrix blk;
//        DenseMatrix subBlk(acc*(dims[i]-2),acc*dims[i]);  subBlk.setZero();
//        for (int j = 0; j < dims[i]-2; j++)
//        {
//            DenseMatrix Ij(acc, acc); Ij.setIdentity();
//            subBlk.block(acc*j,acc*j,acc,3*acc) = kroneckerProduct(d, Ij);
//        }
//        // Complete Kronecker products
//        DenseMatrix tmp1, tmp2;
//        tmp2 = subBlk;

//        for (int j = i+1; j < vars; j++)
//        {
//            tmp1 = tmp2;
//            DenseMatrix Ij(dims[j],dims[j]); Ij.setIdentity();
//            tmp2 = kroneckerProduct(Ij, tmp1);
//        }
//        blk = tmp2;
//        acc *= dims[i];

//        // Insert block into D matrix
//        if (i > 0) insertionRow += numBlkRows[i-1];
//        Dd.block(insertionRow,0,numBlkRows[i],numCols) = blk;
//    }

//    // Convert to sparse matrix and compress
//    D = Dd.sparseView();

    D.makeCompressed();

    return D;
}

void PSplineApproximant::save(const std::string fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void PSplineApproximant::load(const std::string fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

const std::string PSplineApproximant::getDescription() const
{
    std::string description("PSplineApproximant with lambda: ");
    description.append(std::to_string(lambda));
    description.append("\n");
    description.append(bspline.getDescription());

    return description;
}

} // namespace SPLINTER

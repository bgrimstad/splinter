/*
This file is part of the Multivariate Splines library.
Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

#include "generaldefinitions.h"
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseQR"

namespace MultivariateSplines
{

using std::cout;
using std::endl;

template<class lhs, class rhs>
class LinearSolver
{
public:
    bool solve(const lhs &A, const rhs &b, rhs &x) const
    {
        if (!consistentData(A,b))
        {
            cout << "Inconsistent matrix dimensions!" << endl;
            return false;
        }

        bool success = doSolve(A,b,x);

        if (!success || !validSolution(A,b,x))
        {
            cout << "Solver did not converge to acceptable tolerance!" << endl;
            return false;
        }
        return true;
    }

private:
    double tol = 1e-12; // Relative error tolerance

    virtual bool doSolve(const lhs &A, const rhs &b, rhs &x) const = 0;

    bool consistentData(const lhs &A, const rhs &b) const
    {
        if (A.rows() != b.rows())
            return false;
        return true;
    }

    bool validSolution(const lhs &A, const rhs &b, const rhs &x) const
    {
        //return b.isApprox(A*x);
        double err = (A*x - b).norm() / b.norm();
        return (err > tol) ? false : true;
    }
};

class DenseQR : public LinearSolver<DenseMatrix,DenseMatrix>
{
private:
    bool doSolve(const DenseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
    {
        //x = A.colPivHouseholderQr().solve(b);

        // Solve linear system
        Eigen::ColPivHouseholderQR<DenseMatrix> qr(A);
        // Note: qr.info() always returns true
        if (qr.info() == Eigen::Success)
        {
            x = qr.solve(b);
            return true;
        }
        return false;
    }
};

class SparseBiCG : public LinearSolver<SparseMatrix,DenseMatrix>
{
private:
    bool doSolve(const SparseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
    {
        // Init BiCGSTAB solver (requires square matrices)
        Eigen::BiCGSTAB<SparseMatrix> sparseSolver(A);

        if (sparseSolver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparseSolver.solve(b);

            if (sparseSolver.info() == Eigen::Success)
                return true;
        }

        return false;
    }
};

class SparseLU : public LinearSolver<SparseMatrix,DenseMatrix>
{
private:
    bool doSolve(const SparseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
    {
        // Init SparseLU solver (requires square matrices)
        Eigen::SparseLU<SparseMatrix > sparseSolver;
        // Compute the ordering permutation vector from the structural pattern of A
        sparseSolver.analyzePattern(A);
        // Compute the numerical factorization
        sparseSolver.factorize(A);

        if (sparseSolver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparseSolver.solve(b);

            if (sparseSolver.info() == Eigen::Success)
                return true;
        }

        return false;
    }
};

class SparseQR : public LinearSolver<SparseMatrix,DenseMatrix>
{
private:
    bool doSolve(const SparseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
    {
        // Init SparseQR solver (works with rectangular matrices)
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> sparseSolver;
        sparseSolver.analyzePattern(A);
        sparseSolver.factorize(A);

        if (sparseSolver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparseSolver.solve(b);

            if (sparseSolver.info() == Eigen::Success)
                return true;
        }

        return false;
    }
};

} // namespace MultivariateSplines

#endif // LINEARSOLVER_H

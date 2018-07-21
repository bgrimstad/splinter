/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_LINEAR_SOLVER_H
#define SPLINTER_LINEAR_SOLVER_H

#include "definitions.h"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif
#include "Eigen/IterativeLinearSolvers"
#include "Eigen/SparseQR"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif

namespace SPLINTER
{

template<class lhs, class rhs>
class LinearSolver
{
public:
    bool solve(const lhs &A, const rhs &b, rhs &x) const
    {
        if (!consistent_data(A, b))
            throw Exception("LinearSolver::solve: Inconsistent matrix dimensions!");

        bool success = do_solve(A, b, x);

        if (!success)
            throw Exception("LinearSolver::solve: Solver did not converge to acceptable tolerance!");

//        if (!valid_solution(A, b, x))
//            throw Exception("LinearSolver::solve: Invalid solution!");

        return true;
    }

    virtual ~LinearSolver() {}
private:
    double tol = 1e-12; // Relative error tolerance

    virtual bool do_solve(const lhs &A, const rhs &b, rhs &x) const = 0;

    bool consistent_data(const lhs &A, const rhs &b) const
    {
        return A.rows() == b.rows();
    }

    bool valid_solution(const lhs &A, const rhs &b, const rhs &x) const
    {
        //return b.isApprox(A*x);
        double err = (A*x - b).norm() / b.norm();

        return (err <= tol);
    }
};

template<class rhs = DenseVector>
class DenseSVD : public LinearSolver<DenseMatrix, rhs>
{
private:
    bool do_solve(const DenseMatrix &A, const rhs &b, rhs &x) const
    {
        // Solve linear system
        Eigen::JacobiSVD<DenseMatrix> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        x = svd.solve(b);
        return true;
    }
};

template<class rhs = DenseVector>
class DenseQR : public LinearSolver<DenseMatrix, rhs>
{
private:
    bool do_solve(const DenseMatrix &A, const rhs &b, rhs &x) const
    {
        //x = A.colPivHouseholderQr().solve(b);

        // Solve linear system
        Eigen::ColPivHouseholderQR<DenseMatrix> qr(A);

        if (qr.info() == Eigen::Success)
        {
            x = qr.solve(b);

            // Note: qr.info() always returns true
            return qr.info() == Eigen::Success;
        }
        return false;
    }
};

template<class rhs = DenseVector>
class SparseBiCG : public LinearSolver<SparseMatrix, rhs>
{
private:
    bool do_solve(const SparseMatrix &A, const rhs &b, rhs &x) const
    {
        // Init BiCGSTAB solver (requires square matrices)
        Eigen::BiCGSTAB<SparseMatrix> sparse_solver(A);

        if (sparse_solver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparse_solver.solve(b);

            return sparse_solver.info() == Eigen::Success;
        }

        return false;
    }
};

template<class rhs = DenseVector>
class SparseLU : public LinearSolver<SparseMatrix, rhs>
{
private:
    bool do_solve(const SparseMatrix &A, const rhs &b, rhs &x) const
    {
        // Init SparseLU solver (requires square matrices)
        Eigen::SparseLU<SparseMatrix> sparse_solver;
        // Compute the ordering permutation vector from the structural pattern of A
        sparse_solver.analyzePattern(A);
        // Compute the numerical factorization
        sparse_solver.factorize(A);

        if (sparse_solver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparse_solver.solve(b);

            return sparse_solver.info() == Eigen::Success;
        }

        return false;
    }
};

template<class rhs = DenseVector>
class SparseQR : public LinearSolver<SparseMatrix, rhs>
{
private:
    bool do_solve(const SparseMatrix &A, const rhs &b, rhs &x) const
    {
        // Init SparseQR solver (works with rectangular matrices)
        Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> sparse_solver;
        sparse_solver.analyzePattern(A);
        sparse_solver.factorize(A);

        if (sparse_solver.info() == Eigen::Success)
        {
            // Solve LSE
            x = sparse_solver.solve(b);

            return sparse_solver.info() == Eigen::Success;
        }

        return false;
    }
};

} // namespace SPLINTER

#endif // SPLINTER_LINEAR_SOLVER_H

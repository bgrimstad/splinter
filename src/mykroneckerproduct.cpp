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


#include "mykroneckerproduct.h"

namespace MultivariateSplines
{

void myKroneckerProduct(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &AB)
{
    AB.resize(A.rows()*B.rows(), A.cols()*B.cols());

    // Reserve memory for AB

    //AB.reserve(A.nonZeros()*B.nonZeros()); // Does not reserve inner vectors (slow)
    //int innernnz = std::ceil(A.nonZeros()*B.nonZeros()/AB.outerSize());
    //AB.reserve(Eigen::VectorXi::Constant(AB.outerSize(), innernnz)); // Assumes equal distribution of non-zeros (slow)

    // Calculate exact number of non-zeros for each inner vector
    Eigen::VectorXi nnzA = Eigen::VectorXi::Zero(A.outerSize());
    Eigen::VectorXi nnzB = Eigen::VectorXi::Zero(B.outerSize());
    Eigen::VectorXi nnzAB = Eigen::VectorXi::Zero(AB.outerSize());
    //innerNonZeros.setZero();

    for (int jA = 0; jA < A.outerSize(); ++jA)
    {
        int nnz = 0;
        for (SparseMatrix::InnerIterator itA(A,jA); itA; ++itA) nnz++;
        nnzA(jA) = nnz;
    }

    for (int jB = 0; jB < B.outerSize(); ++jB)
    {
        int nnz = 0;
        for (SparseMatrix::InnerIterator itB(B,jB); itB; ++itB) nnz++;
        nnzB(jB) = nnz;
    }

    int innz = 0;
    for (int i = 0; i < nnzA.rows(); ++i)
    {
        for (int j = 0; j < nnzB.rows(); ++j)
        {
            nnzAB(innz) = nnzA(i)*nnzB(j);
            innz++;
        }
    }

    AB.reserve(nnzAB);

    // Non-zero tolerance
    double tolerance = std::numeric_limits<SparseMatrix::Scalar>::epsilon();

    // Compute Kronecker product
    for (int jA = 0; jA < A.outerSize(); ++jA)
    {
        for (SparseMatrix::InnerIterator itA(A,jA); itA; ++itA)
        {
            if (std::abs(itA.value()) > tolerance)
            {
                int jrow = itA.row()*B.rows();
                int jcol = itA.col()*B.cols();

                for (int jB = 0; jB < B.outerSize(); ++jB)
                {
                    for (SparseMatrix::InnerIterator itB(B,jB); itB; ++itB)
                    {
                        if (std::abs(itA.value()*itB.value()) > tolerance)
                        {
                            int row = jrow + itB.row();
                            int col = jcol + itB.col();
                            AB.insert(row,col) = itA.value()*itB.value();
                        }
                    }
                }
            }
        }
    }
    AB.makeCompressed();
    AB.finalize();
}

} // namespace MultivariateSplines

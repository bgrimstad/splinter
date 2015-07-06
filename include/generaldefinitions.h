/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_GENERALDEFINITIONS_H
#define SPLINTER_GENERALDEFINITIONS_H

#ifndef NDEBUG
    #include <iostream>
    #include <iomanip>
#endif // NDEBUG

#ifndef API
# ifdef _MSC_VER
#  define API __declspec(dllexport)
# else
#  define API
# endif
#endif

#include <exception>
#include <stdexcept>

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace SPLINTER
{

// Eigen vectors
typedef Eigen::VectorXd DenseVector;
typedef Eigen::SparseVector<double> SparseVector;

// Eigen matrices
typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::SparseMatrix<double> SparseMatrix; // declares a column-major sparse matrix type of double

// Compare two numbers
template<typename T>
bool assertNear(T x, T y, double tolAbs = 1e-8, double tolRel = 1e-8)
{
    double dx = std::abs(x - y);
    double xAbs = 0.5*(std::abs(x) + std::abs(y));
    double err = std::max(tolAbs, tolRel*xAbs);
    return dx < err;
}

class Exception : public std::exception
{
private:
    std::string __what;

public:

    Exception(const std::string& what)
        : __what(what)
    {
    }

    const char* what() const throw()
    {
        return this->__what.c_str();
    }
};

} // namespace SPLINTER

#endif // SPLINTER_GENERALDEFINITIONS_H
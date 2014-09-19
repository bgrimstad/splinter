/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#ifndef MS_GENERALDEFINITIONS_H
#define MS_GENERALDEFINITIONS_H

#ifndef NDEBUG
    #include <iostream>
    #include <iomanip>
#endif // NDEBUG

#include <exception>
#include <stdexcept>

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace MultivariateSplines
{

// Eigen vectors
typedef Eigen::VectorXd DenseVector;
typedef Eigen::SparseVector<double> SparseVector;

// Eigen matrices
typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::SparseMatrix<double> SparseMatrix; // declares a column-major sparse matrix type of double

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

} // namespace MultivariateSplines

#endif // MS_GENERALDEFINITIONS_H

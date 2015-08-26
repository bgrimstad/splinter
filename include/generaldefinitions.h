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


#ifndef SPLINTER_API
# ifdef _MSC_VER
#  define SPLINTER_API __declspec(dllexport)
# else
#  define SPLINTER_API
# endif
#endif

# ifndef NDEBUG
#  include <iostream>
#  include <iomanip>
# endif // NDEBUG

# include <exception>
# include <stdexcept>

# include <vector>
# include <Eigen/Dense>
# include <Eigen/Sparse>

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

inline std::vector<double> denseVectorToVector(const DenseVector &denseVec)
{
    std::vector<double> vec(denseVec.size());

    for(size_t i = 0; i < denseVec.size(); ++i)
    {
        vec.at(i) = denseVec(i);
    }

    return vec;
}

inline DenseVector vectorToDenseVector(const std::vector<double> &vec)
{
    DenseVector denseVec(vec.size());

    for(size_t i = 0; i < vec.size(); ++i)
    {
        denseVec(i) = vec.at(i);
    }

    return denseVec;
}

inline std::vector<std::vector<double>> denseMatrixToVectorVector(const DenseMatrix &mat)
{
    std::vector<std::vector<double>> vec(mat.rows());

    for(size_t i = 0; i < mat.rows(); ++i)
    {
        for(size_t j = 0; j < mat.cols(); ++j)
        {
            vec.at(i).at(j) = mat(i, j);
        }
    }

    return vec;
}

inline DenseMatrix vectorVectorToDenseMatrix(const std::vector<std::vector<double>> &vec)
{
    DenseMatrix mat(vec.size(), vec.at(0).size());

    for(size_t i = 0; i < vec.size(); ++i)
    {
        for(size_t j = 0; j < vec.at(0).size(); ++j)
        {
            mat(i, j) = vec.at(i).at(j);
        }
    }

    return mat;
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
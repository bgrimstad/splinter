#ifndef GENERALDEFINITIONS_H
#define GENERALDEFINITIONS_H

#include <iostream>
#include <iomanip>

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

} // namespace MultivariateSplines

#endif // GENERALDEFINITIONS_H

#ifndef GENERALDEFINITIONS_H
#define GENERALDEFINITIONS_H

#include <iostream>
#include <iomanip>

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

// Eigen vectors
typedef Eigen::VectorXd DenseVector;
typedef Eigen::SparseVector<double> SparseVector;

// Eigen matrices
typedef Eigen::MatrixXd DenseMatrix;
typedef Eigen::SparseMatrix<double> SparseMatrix; // declares a column-major sparse matrix type of double

// Infinity constant
const double INF = std::numeric_limits<double>::infinity();

enum VariableTypes
{
    CONTINUOUS,
    BINARY,
    INTEGER
};

// Print functions
template<typename T>
void printVector(const std::vector<T> &vec)
{
    for (auto val : vec)
        std::cout << std::setprecision(8) << std::setw(12) << std::left << val << " ";

    std::cout << std::endl;
}

template<typename T>
void printVector(const std::vector< std::vector<T> > &vec)
{
    for (auto val : vec)
        printVector(val);
}

template<typename T>
void printArray(const T a[], int n)
{
    for (int i = 0; i < n; i++)
        std::cout << std::setprecision(4) << std::setw(6) << std::left << a[i] << " ";

    std::cout << std::endl;
}

// Integer test
bool isInteger(double value);

// Randomizer
int randomInteger(int min, int max);

void myKroneckerProduct(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &AB);

#endif // GENERALDEFINITIONS_H

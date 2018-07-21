/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <utils/op_overloads.h>
#include <utils/test_utils.h>
#include <utilities.h>


/*
 * Overload of operators where an Eigen type is one or more of the arguments.
 * Must be defined in the Eigen namespace for Clang to find the overloads.
 */
namespace Eigen
{

bool operator==(const SPLINTER::SparseMatrix &lhs, const SPLINTER::SparseMatrix &rhs)
{
    return SPLINTER::DenseMatrix(lhs) == SPLINTER::DenseMatrix(rhs);
}

bool operator==(const SPLINTER::SparseVector &lhs, const SPLINTER::SparseVector &rhs)
{
    return SPLINTER::DenseVector(lhs) == SPLINTER::DenseVector(rhs);
}

bool operator==(const SPLINTER::DenseVector &denseVec, const std::vector<double> &vec)
{
    return vec == denseVec;
}

bool operator==(const SPLINTER::DenseMatrix &denseMat, const std::vector<std::vector<double>> &vecVec)
{
    return vecVec == denseMat;
}

bool operator==(const std::vector<double> &vec, const SPLINTER::DenseVector &denseVec)
{
    if (vec.size() != (unsigned long) denseVec.size())
        return false;

    for (size_t i = 0; i < vec.size(); ++i)
        if (vec.at(i) != denseVec(i))
            return false;

    return true;
}

bool operator==(const std::vector<std::vector<double>> &vecVec, const SPLINTER::DenseMatrix &denseMat)
{
    size_t matCols = denseMat.cols();
    if(vecVec.size() != (unsigned long) denseMat.rows())
    {
        return false;
    }

    for (size_t i = 0; i < vecVec.size(); ++i)
    {
        if (vecVec.at(i).size() != matCols)
            return false;

        for (size_t j = 0; j < matCols; ++j)
            if (vecVec.at(i).at(j) != denseMat(i, j))
                return false;
    }

    return true;
}

}

namespace SPLINTER
{

/*
 * Comparison operators (==)
 */
bool operator==(const DataTable &lhs, const DataTable &rhs)
{
    // NOTE: requires equal ordering of samples and an equal number of samples
    return  lhs._dim_x == rhs._dim_x
            && lhs._dim_y == rhs._dim_y
            && lhs.samples == rhs.samples
            && lhs.get_num_samples() == rhs.get_num_samples();
}

bool operator==(const DataPoint &lhs, const DataPoint &rhs)
{
    for (unsigned int i = 0; i < lhs.get_dim_x(); i++)
    {
        if (!assert_near(lhs.get_x().at(i), rhs.get_x().at(i)))
            return false;
    }

    for (unsigned int i = 0; i < lhs.get_dim_y(); i++)
    {
        if (!assert_near(lhs.get_y().at(i), rhs.get_y().at(i)))
            return false;
    }

    return true;
}

bool operator==(const BSpline &lhs, const BSpline &rhs)
{
    return
            lhs.dim_x == rhs.dim_x
            && lhs.dim_y == rhs.dim_y
            && lhs.control_points.isApprox(rhs.control_points)
            && lhs.basis == rhs.basis
            && lhs.get_num_basis_functions_per_variable() == rhs.get_num_basis_functions_per_variable()
            && lhs.get_knot_vectors() == rhs.get_knot_vectors()
            && lhs.get_basis_degrees() == rhs.get_basis_degrees()
            && lhs.get_domain_lower_bound() == rhs.get_domain_lower_bound()
            && lhs.get_domain_upper_bound() == rhs.get_domain_upper_bound();
}

bool operator==(const BSplineBasis &lhs, const BSplineBasis &rhs)
{
    return
            lhs.num_variables == rhs.num_variables
            && lhs.bases == rhs.bases;
}

bool operator==(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs)
{
    return
            lhs.degree == rhs.degree
            && lhs.knots == rhs.knots
            && lhs.target_num_basis_functions == rhs.target_num_basis_functions;
}

/*
 * Comparison operators (!=)
 */
bool operator!=(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs)
{
    return !(lhs == rhs);
}

bool operator!=(const DataPoint &lhs, const DataPoint &rhs)
{
	return !(lhs == rhs);
}

/*
 * Output stream operator
 */
std::ostream &operator<<(std::ostream &out, const DataPoint &sample) {
    out << "(";
    bool firstLoop = true;
    for (auto val : sample.get_x()) {
        if (!firstLoop) {
            out << ", ";
        }
        out << val;
        firstLoop = false;
    }
    out << ") = (" << sample.get_y() << ")";

    return out;
}

std::ostream &operator<<(std::ostream &out, const DataTable &table)
{
    out << "DataTable (";
    out << "dim_x: " << table.get_dim_x();
    out << ", dim_y: " << table.get_dim_y();
    out << ", num_samples: " << table.get_num_samples();
    out << ")" << std::endl;
    out << "--------------------------------------------------------" << std::endl;

    for (auto &point : table.get_samples()) {
        for (auto x : point.get_x()) {
            out << x << ", ";
        }
        out << " : ";
        for (auto y : point.get_y()) {
            out << y << ", ";
        }
        out << std::endl;
    }

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &obj)
{
    for (const T &elem : obj)
        out << elem << std::endl;

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &obj)
{
    for (const T &elem : obj)
        out << elem << std::endl;

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::multiset<T> &obj)
{
    for (const T &elem : obj)
        out << elem << std::endl;

    return out;
}

} // namespace SPLINTER

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <operator_overloads.h>
#include <testingutilities.h>

namespace SPLINTER
{

/*
 * Comparison operators
 */
bool operator==(const DataTable &lhs, const DataTable &rhs) {
    return
            lhs.getNumVariables() == rhs.getNumVariables()
            && lhs.getNumSamples() == rhs.getNumSamples()
            && lhs.getSamples() == rhs.getSamples()
            && lhs.getGrid() == rhs.getGrid();
}

bool operator==(const DataSample &lhs, const DataSample &rhs) {
    for(unsigned int i = 0; i < lhs.getDimX(); i++)
    {
        if(!equalsWithinRange(lhs.getX().at(i), rhs.getX().at(i))) {
            return false;
        }
    }

    if(!equalsWithinRange(lhs.getY(), rhs.getY())) {
        return false;
    }

    return true;
}

bool operator==(const BSpline &lhs, const BSpline &rhs)
{
    return
            lhs.getNumVariables() == rhs.getNumVariables()
            && lhs.getNumControlPoints() == rhs.getNumControlPoints()
            && lhs.getNumBasisFunctions() == rhs.getNumBasisFunctions()
            && lhs.getKnotVectors() == rhs.getKnotVectors()
            && lhs.getBasisDegrees() == rhs.getBasisDegrees()
            && lhs.getDomainLowerBound() == rhs.getDomainLowerBound()
            && lhs.getDomainUpperBound() == rhs.getDomainUpperBound();
}

template <class T>
bool operator==(const std::vector<T> &lhs, const std::vector<T> &rhs)
{
    auto lit = lhs.cbegin(), rit = rhs.cbegin();
    for (; lit != lhs.cend() && rit != rhs.cend(); ++lit, ++rit)
    {
        if(*lit != *rit) {
            return false;
        }
    }
    return true;
}

template <class T>
bool operator==(const std::set<T> &lhs, const std::set<T> &rhs)
{
    auto lit = lhs.cbegin(), rit = rhs.cbegin();
    for (; lit != lhs.cend() && rit != rhs.cend(); ++lit, ++rit)
    {
        if(*lit != *rit) {
            return false;
        }
    }
    return true;
}

template <class T>
bool operator==(const std::multiset<T> &lhs, const std::multiset<T> &rhs)
{
    auto lit = lhs.cbegin(), rit = rhs.cbegin();
    for (; lit != lhs.cend() && rit != rhs.cend(); ++lit, ++rit)
    {
        if(*lit != *rit) {
            return false;
        }
    }
    return true;
}

template <class T>
bool operator!=(const T &lhs, const T &rhs) {
    return !(lhs == rhs);
}

/*
 * Output stream operator
 */
std::ostream &operator<<(std::ostream &out, const DataSample &sample) {
    out << "(";
    bool firstLoop = true;
    for(auto val : sample.getX()) {
        out << val;
        if(!firstLoop) {
            out << ", ";
        }
        firstLoop = false;
    }
    out << ") = (" << sample.getY() << ")";

    return out;
}

std::ostream &operator<<(std::ostream &out, const DataTable &table) {
    out << "numVariables: " << table.getNumVariables() << std::endl;
    out << "numSamples: " << table.getNumSamples() << std::endl;
    // out << "samples: " << table.getSamples() << std::endl;
    out << "grid dimensions: ";
    bool firstLoop = true;
    for(const auto &dimension : table.getGrid()) {
        out << dimension.size();
        if(!firstLoop) {
            out << ", ";
        }
        firstLoop = false;
    }

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &obj) {
    for(const T &elem : obj) {
        out << elem << std::endl;
    }

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::set<T> &obj) {
    for(const T &elem : obj) {
        out << elem << std::endl;
    }

    return out;
}

template <class T>
std::ostream &operator<<(std::ostream &out, const std::multiset<T> &obj) {
    for(const T &elem : obj) {
        out << elem << std::endl;
    }

    return out;
}

} // namespace SPLINTER

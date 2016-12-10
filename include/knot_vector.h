/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_KNOT_VECTOR_H
#define SPLINTER_KNOT_VECTOR_H

#include <vector>

namespace SPLINTER
{

class KnotVector
{
public:
    KnotVector(const std::vector<double> &knots, unsigned int degree = 0)
            : knots(std::vector<double>(knots)),
              degree(degree)
    {
        // Test knots here

    }

    std::vector<double>::size_type size() const {
        return knots.size();
    }

    bool empty() const {
        return knots.empty();
    }

    std::vector<double>::const_iterator cbegin() const {
        return knots.cbegin();
    }

    std::vector<double>::const_iterator cend() const {
        return knots.cend();
    }

    const double at(unsigned int i) {
        return knots.at(i);
    }

    const double front() {
        return knots.front();
    }

    const double back() {
        return knots.back();
    }

private:
    std::vector<double> knots;
    const unsigned int degree;

//    friend class Serializer;
//    friend bool operator==(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs);
//    friend bool operator!=(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs);
};

bool operator==(const KnotVector &lhs, const KnotVector &rhs);
bool operator!=(const KnotVector &lhs, const KnotVector &rhs);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_VECTOR_H

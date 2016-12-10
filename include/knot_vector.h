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
#include <algorithm>

namespace SPLINTER
{

class KnotVector
{
public:
    KnotVector(const std::vector<double> &knots)
            : knots(std::vector<double>(knots))
    {
        // Test knots here

    }

    std::vector<double> get_raw() const
    {
        // Return copy of knots
        return knots;
    }

    bool is_regular(unsigned int degree);

    bool is_clamped(unsigned int degree);

    bool is_refinement(const std::vector<double> &refinedKnots);

    bool inside_support(double x) const
    {
        return (knots.front() <= x) && (x <= knots.back());
    }

    unsigned int knot_multiplicity(double tau) const
    {
        return (unsigned int)std::count(knots.begin(), knots.end(), tau);
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

//    friend class Serializer;
};

bool operator==(const KnotVector &lhs, const KnotVector &rhs);
bool operator!=(const KnotVector &lhs, const KnotVector &rhs);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_VECTOR_H

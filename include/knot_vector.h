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
#include "definitions.h"

namespace SPLINTER
{

/**
 * Class representing a knot vector (nondecreasing sequence of number)
 */
class KnotVector
{
public:
    KnotVector()
        : knots(std::vector<double>())
    {
    }

    KnotVector(const std::vector<double> &knots)
        : knots(std::vector<double>(knots))
    {
        // Test knots here
        if (!is_nondecreasing())
            throw Exception("KnotVector::KnotVector: Knot vector is not nondecreasing.");
    }

    std::vector<double> get_values() const {
        // Return copy of knots
        return knots;
    }

    bool is_nondecreasing() const {
        return std::is_sorted(knots.begin(), knots.end());
    }

    bool is_regular(unsigned int degree) const;

    bool is_clamped(unsigned int degree) const;

    bool is_refinement(const KnotVector &refined_knots) const {
        return is_refinement(refined_knots.get_values());
    }

    bool is_refinement(const std::vector<double> &refined_knots) const;

    bool is_supported(double x) const {
        return (knots.front() <= x) && (x <= knots.back());
    }

    unsigned int multiplicity(double tau) const {
        return (unsigned int)std::count(knots.begin(), knots.end(), tau);
    }

    unsigned int index_interval(double x) const;

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

    const double at(unsigned int i) const {
        return knots.at(i);
    }

    const double front() const {
        return knots.front();
    }

    const double back() const {
        return knots.back();
    }

private:
    std::vector<double> knots;

    friend class Serializer;
};

bool operator==(const KnotVector &lhs, const KnotVector &rhs);
bool operator!=(const KnotVector &lhs, const KnotVector &rhs);

} // namespace SPLINTER

#endif // SPLINTER_KNOT_VECTOR_H

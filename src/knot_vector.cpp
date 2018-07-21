/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <knot_vector.h>
#include <algorithm>


namespace SPLINTER
{

bool KnotVector::is_regular(unsigned int degree) const
{
//    // Check size
//    if (size() < 2 * (degree + 1))
//        return false;

    // Check order
    if (!is_nondecreasing())
        return false;

    // Check multiplicity of knots
    for (auto it = cbegin(); it != cend(); ++it)
    {
        if (count(cbegin(), cend(), *it) > degree + 1)
            return false;
    }

    return true;
}

bool KnotVector::is_refinement(const std::vector<double> &refined_knots) const
{
    // Check size
    if (refined_knots.size() < knots.size())
        return false;

    // Check that each element in knots occurs at least as many times in refined_knots
    for (auto it = knots.cbegin() ; it != knots.cend(); ++it)
    {
        auto m_tau = count(knots.begin(), knots.end(), *it);
        auto m_t = count(refined_knots.begin(), refined_knots.end(), *it);
        if (m_t < m_tau) return false;
    }

    // Check that range is not changed
    if (knots.front() != refined_knots.front() || knots.back() != refined_knots.back())
        return false;

    return true;
}

/**
 * Check if end knots have multiplicity degree + 1
 * @param degree
 * @return bool
 */
bool KnotVector::is_clamped(unsigned int degree) const
{
    // Check multiplicity of first knot
    if (std::count(knots.begin(), knots.begin() + degree + 1, knots.front()) != degree + 1)
        return false;

    // Check multiplicity of last knot
    if (std::count(knots.end() - degree - 1, knots.end(), knots.back()) != degree + 1)
        return false;

    return true;
}

/**
 * Finds index i such that knots.at(i) <= x < knots.at(i+1).
 */
unsigned int KnotVector::index_interval(double x) const
{
    if (!is_supported(x))
        throw Exception("KnotVector::index_interval: x outside knot support!");

    // Find first knot that is larger than x
    auto it = std::upper_bound(cbegin(), cend(), x);

    // Compute index
    auto index = (it - knots.cbegin()) - 1;

    if (index < 0)
        throw Exception("KnotVector::index_interval: computed negative index!");

    return (unsigned int)index;
}

bool operator==(const KnotVector &lhs, const KnotVector &rhs) {
    auto is_equal = std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
    return is_equal && lhs.size() == rhs.size();
}

bool operator!=(const KnotVector &lhs, const KnotVector &rhs) {
    return !(lhs == rhs);
}

} // namespace SPLINTER

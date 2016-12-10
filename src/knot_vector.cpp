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

bool operator==(const KnotVector &lhs, const KnotVector &rhs) {
    auto is_equal = std::equal(lhs.cbegin(), lhs.cend(), rhs.cbegin());
    return is_equal && lhs.size() == rhs.size();
}

bool operator!=(const KnotVector &lhs, const KnotVector &rhs) {
    return !(lhs == rhs);
}

} // namespace SPLINTER

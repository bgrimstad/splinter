/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <bspline.h>
#include <bspline_utils.h>

using namespace SPLINTER;

#define COMMON_TAGS "[unit][bsplineutils]"
#define COMMON_TEXT " unit test"

TEST_CASE("computeSecondOrderDifferenceMatrix 1-D" COMMON_TEXT, COMMON_TAGS)
{
    std::vector<double> knots = {1, 1, 1, 1, 2, 3, 4, 5, 5, 5, 5};
    std::vector<std::vector<double>> knot_vectors = {knots};
    std::vector<unsigned int> degrees = {3};

    BSpline bs(degrees, knot_vectors);

    auto D = compute_second_order_finite_difference_matrix(bs);

    unsigned int n = bs.get_num_basis_functions();
    SparseMatrix D_expected(n-2, n);

    for (unsigned int i = 0; i < n - 2; ++i)
    {
        D_expected.insert(i, i) = 1;
        D_expected.insert(i, i+1) = -2;
        D_expected.insert(i, i+2) = 1;
    }

    REQUIRE(D.isApprox(D_expected));
}

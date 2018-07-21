/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINETESTINGUTILITIES_H
#define SPLINTER_BSPLINETESTINGUTILITIES_H

#include <bspline_builders.h>
#include <data_table.h>

namespace SPLINTER
{

DataTable sample_test_function();

/*
 * Test knot insertion
 */
bool test_knot_insertion();

/*
 * Methods for B-spline domain reduction testing
 */
bool domain_reduction_test(BSpline &bs, const BSpline &bs_orig);
bool run_recursive_domain_reduction_test();
bool domain_reduction_test1();

} // namespace SPLINTER

#endif //SPLINTER_BSPLINETESTINGUTILITIES_H

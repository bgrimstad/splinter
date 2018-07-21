/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <utilities.h>
#include "json_parser.h"
#include "utils/bspline_collection.h"

using namespace SPLINTER;

#define COMMON_TAGS "[unit][json]"

TEST_CASE("BSpline can be saved and loaded from json", COMMON_TAGS)
{
    std::string filename = "bspline.json";

    auto col = get_bspline_collection();

    for (const auto &bspline : col) {
        bspline.to_json(filename);
        BSpline loaded_bspline = BSpline::from_json(filename);
        REQUIRE(bspline == loaded_bspline);
    }

    remove(filename.c_str());
}

TEST_CASE("DataTable can be saved and loaded from json", COMMON_TAGS)
{
    std::string filename = "datatable.json";

    auto datatable = DataTable();

    for (auto x1 : linspace(1, 20, 20)) {
        for (auto x2 : linspace(1, 20, 20)) {
            std::vector<double> x = {x1, x2};
            std::vector<double> y = {x1+x2, 1.00001*x1*x2};
            datatable.add_sample(x, y);
        }
    }

    datatable.to_json(filename);
    auto loaded_datatable = DataTable::from_json(filename);
    REQUIRE(datatable == loaded_datatable);

    remove(filename.c_str());
}


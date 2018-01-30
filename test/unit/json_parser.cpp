/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include "json_parser.h"
#include "utils/bspline_collection.h"

using namespace SPLINTER;


#define COMMON_TAGS "[unit][json]"


TEST_CASE("BSpline can be saved and loaded from json", COMMON_TAGS)
{
    std::string filename = "bspline.json";

    auto col = get_bspline_collection();

    for (const auto &bspline : col) {
        save_to_json(bspline, filename);
        BSpline loadedBspline = load_from_json(filename);
        REQUIRE(bspline == loadedBspline);
    }

    remove(filename.c_str());
}

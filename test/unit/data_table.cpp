/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <iostream>
#include <utils/test_utils.h>
#include <utils/test_function_utils.h>
#include <data_table.h>
#include <utilities.h>
#include <utils/data_table_collection.h>

using namespace SPLINTER;

#define COMMON_TAGS "[general][datatable]"

TEST_CASE("DataTable save and load data tables", COMMON_TAGS) {
    auto collection = get_data_table_collection();

    for (auto &table : collection)
    {
        auto dim_x = table.get_dim_x();
        auto dim_y = table.get_dim_y();
        std::string filename = std::string("data_table_")
                               + std::to_string(dim_x)
                               + std::string("_")
                               + std::to_string(dim_y)
                               + std::string(".json");

        table.to_json(filename);

        // Load table and compare
        auto new_table = DataTable::from_json(filename);
        CHECK(new_table == table);

        // Delete saved file
        remove(filename.c_str());
    }
}

TEST_CASE("DataTable load data tables", COMMON_TAGS) {

    // Files to load must correspond to the collection of DataTables
    auto collection = get_data_table_collection();

    std::string dir = "test-resources/data_tables/";
    std::string dir_old_format = "test-resources/data_tables_old_format/";

    for (auto &table : collection)
    {
        auto dim_x = table.get_dim_x();
        auto dim_y = table.get_dim_y();

        std::string filename = std::string("data_table_")
                               + std::to_string(dim_x)
                               + std::string("_")
                               + std::to_string(dim_y)
                               + std::string(".json");

        // Old format
        auto loaded_table_old_format = DataTable::from_json(dir_old_format + filename);
        CHECK((loaded_table_old_format.get_dim_x() == dim_x && loaded_table_old_format.get_dim_y() == dim_y));

        // New format
        auto loaded_table = DataTable::from_json(dir + filename);
        CHECK((loaded_table.get_dim_x() == dim_x && loaded_table.get_dim_y() == dim_y));
    }
}

TEST_CASE("DataTable add samples", COMMON_TAGS) {
    DataTable table, table2, table3, table4, table5;

    DataPoint x1({1, 2}, 3);
    DataPoint x2({1.1, 2}, 3);
    DataPoint x3({1.1, 2.2}, 3);

    table.add_sample(x1);
    table.add_sample(x2);
    table.add_sample(x3);

    table2.add_sample({ {1, 2}, 3 });
    table2.add_sample({ {1.1, 2}, 3 });
    table2.add_sample({ {1.1, 2.2}, 3 });

    // Different order
    table3.add_sample(x3);
    table3.add_sample(x2);
    table3.add_sample(x1);

    CHECK(table == table2);
    CHECK(!(table == table3));

    // Using vectors
    std::vector<double> x1_vec = {0.1};
    std::vector<double> x2_vec = {0.2};
    std::vector<double> x3_vec = {0.3};

    std::vector<double> y1_vec = {0.1, 0.2};
    std::vector<double> y2_vec = {0.2, 0.3};
    std::vector<double> y3_vec = {0.3, 0.4};

    table4.add_sample(x1_vec.at(0), y1_vec);
    table4.add_sample(x2_vec.at(0), y2_vec);
    table4.add_sample(x3_vec.at(0), y3_vec);

    table5.add_sample(x1_vec, y1_vec);
    table5.add_sample(x2_vec, y2_vec);
    table5.add_sample(x3_vec, y3_vec);

    CHECK(table4 == table5);
}

TEST_CASE("DataTable initializer_list behaviour", COMMON_TAGS)
{
    std::vector<DataPoint> list{ { 1, 1 }, { 1, 1 } };
    DataTable table, table_ref;
    table.add_sample({ { 1,1 }, { 1,1 } });
    for (auto& sample : list)
    {
        table_ref.add_sample(sample);
    }
    REQUIRE(table.get_samples() == table_ref.get_samples());
}


TEST_CASE("DataTable copy constructor", COMMON_TAGS) {
    int dim = 2;
    auto func = getTestFunction(dim, 2);

    std::vector<double> start = {0, 0};
    std::vector<double> end = {9.0, 4.0};
    std::vector<unsigned int> num_points = {1000, 1000};
    unsigned int num_samples = 1000*1000;

    // range1 is from 0.0,0.0 to 9.0,4.0, 50 points total
    auto points = multi_linspace(start, end, num_points);

    DataTable table;

    for(auto &p : points) {
        table.add_sample(p, func->eval(p));
    }

    CHECK(table.get_num_samples() == num_samples);

    auto table2 = DataTable(table);
    CHECK(table == table2);
}

TEST_CASE("DataTable is_grid_complete", COMMON_TAGS) {

    auto x1 = linspace(0, 1.0, 10);
    auto x2 = linspace(0, 1.0, 10);

    DataTable table, table2;

    int counter = 0;
    for(auto &x1_i : x1) {
        for (auto &x2_i : x2) {
            table.add_sample({x1_i, x2_i}, 0.);
            // Do not add all samples to table2
            if (counter != 1) {
                table2.add_sample({x1_i, x2_i}, 0.);
            }

            ++counter;
        }
    }

    CHECK(table.is_grid_complete());
    CHECK(!table2.is_grid_complete());
}

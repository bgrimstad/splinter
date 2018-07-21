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
#include <data_table.h>
#include <data_table2.h>
#include <utils/data_table_collection.h>


using namespace SPLINTER;


#define COMMON_TAGS "[general][datatable]"


DataTable create_datatable_old_format() {
    int dim = 2;
    auto func = getTestFunction(dim, 2);

    std::vector<double> start = {0, 0};
    std::vector<double> end = {9.0, 4.0};
    std::vector<unsigned int> num_points = {10, 5};

    // Range [0.0, 9.0] x [0.0, 4.0], 50 points total
    auto points = multi_linspace(start, end, num_points);

    DataTable table;

    for(auto &p : points) {
        table.add_sample(p, func->eval(p));
    }

    return table;
}

DataTable2 copy_datatable(DataTable table) {

    DataTable2 new_table;

    for (auto &point : table.get_samples()) {
        new_table.add_sample(point);
    }

    return new_table;
}

bool compare_datatables(DataTable table1, DataTable2 table2) {

    if (table1.get_dim_x() != table2.get_dim_x()) {
        return false;
    }

    if (table1.get_dim_y() != table2.get_dim_y()) {
        return false;
    }

    if (table1.get_num_samples() != table2.get_num_samples()) {
        return false;
    }

    auto it1 = table1.cbegin();
    auto it2 = table2.cbegin();

    while (it1 != table1.cend()) {
        auto p1 = *it1;
        auto p2 = *it2;

        if (p1 != p2) {
            return false;
        }

        ++it1;
        ++it2;
    }

    return true;
}

TEST_CASE("DataTable save data tables", COMMON_TAGS) {
    auto collection = get_data_table_collection();

    for (auto &table : collection)
    {
        auto dim_x = table.get_dim_x();
        auto dim_y = table.get_dim_y();
        std::string filename = std::string("data_table_")
                               + std::to_string(dim_x)
                               + std::string("_")
                               + std::to_string(dim_y)
                               + std::string("_old_format")
                               + std::string(".json");

        table.to_json(filename);

        // Load table and compare
        auto new_table = DataTable2::from_json(filename);
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
        auto loaded_table_old_format = DataTable2::from_json(dir_old_format + filename);
        CHECK((loaded_table_old_format.get_dim_x() == dim_x && loaded_table_old_format.get_dim_y() == dim_y));

        // New format
        auto loaded_table = DataTable2::from_json(dir + filename);
        CHECK((loaded_table.get_dim_x() == dim_x && loaded_table.get_dim_y() == dim_y));
    }
}

TEST_CASE("DataTable2 add samples", COMMON_TAGS) {
    DataTable2 table, table2, table3;

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
}

TEST_CASE("DataTable initializer_list behaviour", COMMON_TAGS)
{
    std::vector<DataPoint> list{ { 1, 1 }, { 1, 1 } };
    DataTable2 table, table_ref;
    table.add_sample({ { 1,1 }, { 1,1 } });
    for (auto& sample : list)
    {
        table_ref.add_sample(sample);
    }
    REQUIRE(table.get_samples() == table_ref.get_samples());
}


TEST_CASE("DataTable2 addition", COMMON_TAGS) {
    int dim = 2;
    auto func = getTestFunction(dim, 2);

    std::vector<double> start = {0, 0};
    std::vector<double> end = {9.0, 4.0};
    std::vector<unsigned int> num_points = {10, 5};

    // range1 is from 0.0,0.0 to 9.0,4.0, 50 points total
    auto points = multi_linspace(start, end, num_points);

    DataTable table, table2;

    for(auto &p : points) {
        table.add_sample(p, func->eval(p));
        table2.add_sample(p, func->eval(p));
    }

    for (auto &sample : table.get_samples())
    {
        for (auto x : sample.get_x()) {
            std::cout << x << ", ";
        }
        std::cout << " : ";
        for (auto y : sample.get_y()) {
            std::cout << y << ", ";
        }
        std::cout << std::endl;
    }

    DataPoint x1({1, 2}, 3);
    DataPoint x2({1.1, 2}, 3);
    DataPoint x3({1.1, 2.2}, 3);

    std::vector<DataPoint> table3, table4;
    table3.emplace_back(x1);
    table3.emplace_back(x2);
    table3.emplace_back(x3);
    table4.emplace_back(x1);
    table4.emplace_back(x2);
//    table4.emplace_back(x3);

    CHECK(table3 == table4);

}

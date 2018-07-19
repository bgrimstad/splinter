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

using namespace SPLINTER;


#define COMMON_TAGS "[general][datatable2]"


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

TEST_CASE("DataTable2 add points", COMMON_TAGS) {
    DataTable2 table, table2;

    DataPoint x1({1, 2}, 3);
    DataPoint x2({1.1, 2}, 3);
    DataPoint x3({1.1, 2.2}, 3);

    table.add_sample(x1);
    table.add_sample(x2);
    table.add_sample(x3);

    table2.add_sample(x1);
    table2.add_sample(x2);
    table2.add_sample(x3);

//    table2.add_sample({1, 2}, 3);
//    table2.add_sample({1, 2}, 3);
//    table2.add_sample({1, 2}, 3);

    CHECK(table == table2);

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

TEST_CASE("DataTable2 save and load", COMMON_TAGS) {
    int dim = 2;
    auto func = getTestFunction(dim, 2);

    std::vector<double> start = {0, 0};
    std::vector<double> end = {9.0, 4.0};
    std::vector<unsigned int> num_points = {10, 5};

    // range1 is from 0.0,0.0 to 9.0,4.0, 50 points total
    auto points = multi_linspace(start, end, num_points);

    DataTable table;
    DataTable2 table2;

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

    std::string filename = "data.json";

    table.to_json(filename);
    auto table5 = DataTable2::from_json(filename);
    std::cout << "First table" << std::endl;
    std::cout << table << std::endl;
    std::cout << "Second table" << std::endl;
    std::cout << table5 << std::endl;

    remove(filename.c_str());

    CHECK(compare_datatables(table, table5));
}

TEST_CASE("DataTable2 set operations", COMMON_TAGS) {
    int dim = 2;
    auto func = getTestFunction(dim, 2);

    auto start = std::vector<double>(dim);
    auto end = std::vector<double>(dim);
    auto points = std::vector<unsigned int>(dim);

    /* Grid:
     *   0 1 2 3 4 5 6 7 8 9
     * 0 x x x x x x x x x x
     * 1 x x x x x x x x x x
     * 2 x x x x x x x x x x
     * 3 x x x x x x x x x x
     * 4 x x x x x x x x x x
     * 5 o o o o o o o o o o
     * 6 o o o o o o o o o o
     * 7 o o o o o o o o o o
     * 8 o o o o o o o o o o
     * 9 o o o o o o o o o o
     *
     * x: samples in table1
     * o: samples in table2
     */

    start.at(0) = 0.0;
    start.at(1) = 0.0;
    end.at(0) = 9.0;
    end.at(1) = 4.0;
    points.at(0) = 10;
    points.at(1) = 5;
    // range1 is from 0.0,0.0 to 9.0,4.0, 50 points total
    auto range1 = multi_linspace(start, end, points);
    auto table1 = sample(func, range1);

    start.at(0) = 0.0;
    start.at(1) = 5.0;
    end.at(0) = 9.0;
    end.at(1) = 9.0;
    points.at(0) = 10;
    points.at(1) = 5;
    // range2 is from 0.0,5.0 to 9.0,9.0, 50 points total
    auto range2 = multi_linspace(start, end, points);
    auto table2 = sample(func, range2);

    start.at(0) = 0.0;
    start.at(1) = 0.0;
    end.at(0) = 9.0;
    end.at(1) = 9.0;
    points.at(0) = 10;
    points.at(1) = 10;
    // range3 is from 0.0,0.0 to 9.0,9.0, 100 points total
    auto range3 = multi_linspace(start, end, points);
    auto table3 = sample(func, range3);

//    // Summation
//    CHECK(table1 + table2 == table3);
//    CHECK(table2 + table1 == table3);
//
//    // Subtraction
//    CHECK(table3 - table2 == table1);
//    CHECK(table3 - table1 == table2);
//
//    // A table subtracted from itself should yield a table with no samples
//    DataTable zeroSampleTable;
//    CHECK(table3 - table3 == zeroSampleTable);
//    CHECK(table3 - table2 - table1 == zeroSampleTable);
//
//    // A table added to itself should be equal to itself (disregarding allowing duplicates)
//    auto table4 = table3 + table3;
//    CHECK(table4 == table3);
}

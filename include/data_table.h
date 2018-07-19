/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_DATA_TABLE_H
#define SPLINTER_DATA_TABLE_H

#include <set>
#include "data_point.h"
#include <ostream>


namespace SPLINTER
{

/*
 * DataTable is a class for storing multidimensional data samples (x, y).
 * The samples are stored in a continuously sorted table.
 */
class SPLINTER_API DataTable
{
public:
    DataTable();
    DataTable(bool allow_duplicates);
    DataTable(bool allow_duplicates, bool allow_incomplete_grid);

    /*
     * Functions for adding a sample (x, y)
     */
    void add_sample(const DataPoint &sample);
    void add_sample(double x, double y);
    void add_sample(const std::vector<double> &x, double y);
    void add_sample(double x, const std::vector<double> &y);
    void add_sample(const std::vector<double> &x, const std::vector<double> &y);
    void add_sample(std::initializer_list<DataPoint> samples);

    /*
     * Getters
     */
    std::multiset<DataPoint>::const_iterator cbegin() const;
    std::multiset<DataPoint>::const_iterator cend() const;

    unsigned int get_dim_x() const {
        return _dim_x;
    }

    unsigned int get_dim_y() const {
        return _dim_y;
    }

    unsigned int get_num_samples() const {
        return (unsigned int) samples.size();
    }

    const std::multiset<DataPoint>& get_samples() const {
        return samples;
    }

    std::vector<std::set<double>> get_grid() const {
        return grid;
    }

    std::vector<std::vector<double>> get_table_x() const;

    std::vector<std::vector<double>> get_table_y() const;

    bool is_grid_complete() const;

    /**
     * Save and load
     */
    void to_json(const std::string &filename) const;

    static DataTable from_json(const std::string &filename);

private:
    bool _allow_duplicates;
    bool _allow_incomplete_grid;
    unsigned int _num_duplicates;
    unsigned int _dim_x;
    unsigned int _dim_y;

    std::multiset<DataPoint> samples;
    std::vector<std::set<double>> grid;

    // Initialise grid to be a std::vector of _dim_x sized std::sets
    void init_data_structures();

    unsigned int get_num_samples_required() const;

    void record_grid_point(const DataPoint &sample);

    // Used by functions that require the grid to be complete before they start their operation
    // This function prints a message and exits the program if the grid is not complete.
    void grid_complete_guard() const;

    friend bool operator==(const DataTable &lhs, const DataTable &rhs);
    friend void datatable_to_json(const DataTable &data, const std::string &filename);
};

DataTable operator+(const DataTable &lhs, const DataTable &rhs);
DataTable operator-(const DataTable &lhs, const DataTable &rhs);

} // namespace SPLINTER

#endif // SPLINTER_DATA_TABLE_H

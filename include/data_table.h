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

#include "data_point.h"
#include <ostream>
#include <json_parser.h>


namespace SPLINTER
{

/*
 * DataTable is a class for storing multidimensional data samples (x, y).
 * The samples are stored in a vector with properties:
 * - Elements can be stored in any order (no sorted order)
 * - Duplicates are allowed (DataPoint x-values determines uniqueness)
 */
class SPLINTER_API DataTable
{
public:
    DataTable() : _dim_x(1),  _dim_y(1) {}

    /*
     * Functions for adding a sample (x, y)
     */
    void add_sample(const DataPoint &sample);

    void add_sample(std::initializer_list<DataPoint> samples);

    template <class Tx, class Ty>
    void add_sample(Tx x, Ty y) {
        return add_sample(DataPoint(x, y));
    }

    template <class Tx, class Ty>
    void add_sample(std::initializer_list<Tx> x, Ty y) {
        return add_sample(DataPoint(x, y));
    }

    template <class Tx, class Ty>
    void add_sample(Tx x, std::initializer_list<Ty> y) {
        return add_sample(DataPoint(x, y));
    }

    template <class Tx, class Ty>
    void add_sample(std::initializer_list<Tx> x, std::initializer_list<Ty> y) {
        return add_sample(DataPoint(x, y));
    }

    /*
     * Getters
     */
    std::vector<DataPoint>::const_iterator cbegin() const {
        return samples.cbegin();
    }

    std::vector<DataPoint>::const_iterator cend() const {
        return samples.cend();
    }

    unsigned int get_dim_x() const {
        return _dim_x;
    }

    unsigned int get_dim_y() const {
        return _dim_y;
    }

    unsigned int get_num_samples() const {
        return (unsigned int) samples.size();
    }

    const std::vector<DataPoint>& get_samples() const {
        return samples;
    }

    std::vector<std::vector<double>> get_table_x() const;

    std::vector<std::vector<double>> get_table_y() const;

    /**
     * Save and load
     */
    void to_json(const std::string &filename) const {
        SPLINTER::datatable_to_json(*this, filename);
    }

    static DataTable from_json(const std::string &filename) {
        return SPLINTER::datatable_from_json(filename);
    }

    /**
     * Utilities
     */
    bool is_grid_complete() const;

private:
    unsigned int _dim_x;
    unsigned int _dim_y;
    std::vector<DataPoint> samples;

    friend bool operator==(const DataTable &lhs, const DataTable &rhs);
//    friend void datatable_to_json(const DataTable &data, const std::string &filename);
};

} // namespace SPLINTER

#endif // SPLINTER_DATA_TABLE_H

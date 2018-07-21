/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_DATA_POINT_H
#define SPLINTER_DATA_POINT_H

#include "definitions.h"

namespace SPLINTER
{

/*
 * DataPoint is a class representing a data point (x, y), where y is the value obtained by sampling at a point x.
 * Note that x and y may be a scalar or vector.
 */
class DataPoint
{
public:
    DataPoint(double x, double y);
    DataPoint(const std::vector<double> &x, double y);
    DataPoint(double x, const std::vector<double> &y);
    DataPoint(const std::vector<double> &x, const std::vector<double> &y);

    /**
     * Getters
     */
    std::vector<double> get_x() const {
        return x;
    }

    std::vector<double> get_y() const {
        return y;
    }

    unsigned int get_dim_x() const {
        return (unsigned int) x.size();
    }

    unsigned int get_dim_y() const {
        return (unsigned int) y.size();
    }

    bool operator<(const DataPoint &rhs) const; // Returns false if the two are equal

private:
    DataPoint() {};

    std::vector<double> x;
    std::vector<double> y;
    void set_data(const std::vector<double> &x, const std::vector<double> &y);
};

} // namespace SPLINTER

#endif // SPLINTER_DATA_POINT_H

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "data_point.h"

namespace SPLINTER
{

DataPoint::DataPoint()
{
}

DataPoint::DataPoint(double x, double y)
{
    setData(std::vector<double>(1, x),
            std::vector<double>(1, y));
}

DataPoint::DataPoint(const std::vector<double> &x, double y)
{
    setData(x, std::vector<double>(1, y));
}

DataPoint::DataPoint(double x, const std::vector<double> &y)
{
    setData(std::vector<double>(1, x), y);
}

DataPoint::DataPoint(const std::vector<double> &x, const std::vector<double> &y)
{
    setData(x, y);
}

void DataPoint::setData(const std::vector<double> &x, const std::vector<double> &y)
{
    this->x = x;
    this->y = y;
}

bool DataPoint::operator<(const DataPoint &rhs) const
{
    if (this->get_dim_x() != rhs.get_dim_x())
        throw Exception("DataPoint::operator<: Cannot compare data points of different dimensions");

    for (unsigned int i = 0; i < this->get_dim_x(); i++)
    {
        if (x.at(i) < rhs.getX().at(i))
            return true;
        else if (x.at(i) > rhs.getX().at(i))
            return false;
    }

    return false;
}

} // namespace SPLINTER

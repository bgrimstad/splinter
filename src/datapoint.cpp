/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "datapoint.h"

namespace SPLINTER
{

DataPoint::DataPoint()
{
}

DataPoint::DataPoint(double x, double y)
{
    setData(std::vector<double>(1, x), y);
}

DataPoint::DataPoint(std::vector<double> x, double y)
{
    setData(x, y);
}

DataPoint::DataPoint(DenseVector x, double y)
{
    std::vector<double> newX;

    for (int i = 0; i < x.size(); i++)
    {
        newX.push_back(x(i));
    }

    setData(newX, y);
}

void DataPoint::setData(const std::vector<double> &x, double y)
{
    this->x = x;
    this->y = y;
}

bool DataPoint::operator<(const DataPoint &rhs) const
{
    if (this->getDimX() != rhs.getDimX())
        throw Exception("DataPoint::operator<: Cannot compare data points of different dimensions");

    for (unsigned int i = 0; i < this->getDimX(); i++)
    {
        if (x.at(i) < rhs.getX().at(i))
            return true;
        else if (x.at(i) > rhs.getX().at(i))
            return false;
    }

    return false;
}

} // namespace SPLINTER

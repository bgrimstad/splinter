/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "datasample.h"

namespace MultivariateSplines
{

DataSample::DataSample(double x, double y)
{
    setData(std::vector<double>(1, x), y);
}

DataSample::DataSample(std::vector<double> x, double y)
{
    setData(x, y);
}

DataSample::DataSample(DenseVector x, double y)
{
    std::vector<double> newX;

    for (int i = 0; i < x.size(); i++)
    {
        newX.push_back(x(i));
    }

    setData(newX, y);
}

void DataSample::setData(const std::vector<double> &x, double y)
{
    this->x = x;
    this->y = y;
}

bool DataSample::operator<(const DataSample &rhs) const
{
    assert(this->getDimX() == rhs.getDimX());

    for(unsigned int i = 0; i < this->getDimX(); i++)
    {
        if(x.at(i) < rhs.getX().at(i))
            return true;
        else if(x.at(i) > rhs.getX().at(i))
            return false;
    }

    return false;
}

std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample)
{
    outputStream << "Sample: (";

    bool firstLoop = true;
    for(auto &coordinate : sample.getX())
    {
        if(!firstLoop)
            outputStream << ", ";

        outputStream << coordinate;
        firstLoop = false;
    }

    outputStream << ") = (" << sample.getY() << ")";

    return outputStream;
}

} // namespace MultivariateSplines

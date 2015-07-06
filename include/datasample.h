/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_DATASAMPLE_H
#define SPLINTER_DATASAMPLE_H

#include "generaldefinitions.h"

namespace SPLINTER
{

/* Class representing a data sample (x,y)
 * where y is the value obtained by sampling
 * at a point x.
*/
class DataSample
{
public:
    DataSample(double x, double y);
    DataSample(std::vector<double> x, double y);
    DataSample(DenseVector x, double y);

    bool operator<(const DataSample &rhs) const; // Returns false if the two are equal
    friend std::ostream &operator<<(std::ostream &outputStream, const DataSample &sample);

    std::vector<double> getX() const { return x; }
    double getY() const { return y; }
    unsigned int getDimX() const { return x.size(); }

private:
    DataSample();

    std::vector<double> x;
    double y;
    void setData(const std::vector<double> &x, double y);

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_DATASAMPLE_H
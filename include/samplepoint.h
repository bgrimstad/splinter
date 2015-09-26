/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_SAMPLEPOINT_H
#define SPLINTER_SAMPLEPOINT_H

#include "definitions.h"

namespace SPLINTER
{

/* Class representing a sample point (x,y), where y is the value
 * obtained by sampling at a point x.
*/
class SamplePoint
{
public:
    SamplePoint(double x, double y);
    SamplePoint(std::vector<double> x, double y);
    SamplePoint(DenseVector x, double y);

    bool operator<(const SamplePoint &rhs) const; // Returns false if the two are equal

    std::vector<double> getX() const { return x; }
    double getY() const { return y; }
    unsigned int getDimX() const { return x.size(); }

private:
    SamplePoint();

    std::vector<double> x;
    double y;
    void setData(const std::vector<double> &x, double y);

    friend class Serializer;
};

} // namespace SPLINTER

#endif // SPLINTER_SAMPLEPOINT_H

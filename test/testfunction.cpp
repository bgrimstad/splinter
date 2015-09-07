/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <testfunction.h>
#include <cmath> // std::ceil

namespace SPLINTER
{

TestFunction::TestFunction(std::function<double (const std::vector<double> &)> f,
                           size_t numVariables,
                           std::string functionString)
        : f(f),
          Function(numVariables),
          functionString(functionString),
          constDegree(false),
          constDegreeVal(std::vector<double>(numVariables, 0.0))
{
}

TestFunction::TestFunction(std::function<double (const std::vector<double> &)> f, size_t numVariables,
                           std::string functionString,  std::vector<double> &degrees)
        : f(f),
          Function(numVariables),
          functionString(functionString),
          constDegree(true),
          constDegreeVal(degrees)
{
}

TestFunction::~TestFunction()
{
}

double TestFunction::eval(const std::vector<double> &x) const
{
    return f(x);
}

std::vector<unsigned int> TestFunction::getConstDegreeInt() const
{
    auto intDegrees = std::vector<unsigned int>(numVariables, 0);

    for (size_t i = 0; i < constDegreeVal.size(); ++i)
    {
        intDegrees.at(i) = (unsigned int) std::ceil(constDegreeVal.at(i));
    }

    return intDegrees;
}

double TestFunction::getMaxDegree() const
{
    double maxDegree = 0.0;
    for (auto deg : constDegreeVal)
    {
        if (deg > maxDegree)
        {
            maxDegree = deg;
        }
    }
    return maxDegree;
}

} // namespace SPLINTER

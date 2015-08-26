/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <testfunction.h>

namespace SPLINTER
{

TestFunction::TestFunction(std::function<double (const std::vector<double> &)> f,
                           size_t numVariables,
                           std::string functionString)
        : f(f),
          Function(numVariables),
          functionString(functionString),
          constDegree(false),
          constDegreeVal(0.0)
{
}

TestFunction::TestFunction(std::function<double (const std::vector<double> &)> f, size_t numVariables,
                           std::string functionString, double constDegreeVal)
        : f(f),
          Function(numVariables),
          functionString(functionString),
          constDegree(true),
          constDegreeVal(constDegreeVal)
{
}

TestFunction::~TestFunction()
{
}

double TestFunction::eval(const std::vector<double> &x) const
{
    return f(x);
}

} // namespace SPLINTER
/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_APPROXIMANT_H
#define SPLINTER_APPROXIMANT_H

#include "definitions.h"
#include "function.h"

namespace SPLINTER
{

/**
 * Interface for approximants on the form:
 * f(x) = c(1)*b1(x) + c(2)*b2(x) + ... + c(n)*bn(x) = c^T * b(x),
 * where c are coefficients and b is a vector of basis functions.
 */
class SPLINTER_API Approximant : public Function
{
public:
    Approximant(unsigned int numVariables)
        : Function(numVariables) {}

    virtual ~Approximant() {}

    /**
     * Serialize and save approximant to fileName
     * Throws if file could not be opened
     */
    virtual void save(const std::string fileName) const = 0;

    /**
     * Deserialize and load approximant from fileName
     * Throws if file could not be opened or if the file format is wrong
     */
    virtual void load(const std::string fileName) = 0;

    /*
     * Functions for appraising absolute approximation error
     */
    //double absDist2(DataTable realValues);
    //double absDistInf(DataTable realValues);

    /*
     * Functions for appraising relative approximation error
     */
    //double relDist2(DataTable realValues);
    //double relDistInf(DataTable realValues);

    /**
     * Returns a name describing this Approximant.
     * Will typically include name of the class and degree
     */
    virtual const std::string getDescription() const = 0;
};

} // namespace SPLINTER

#endif // SPLINTER_APPROXIMANT_H

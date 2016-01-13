/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BUILDER_H
#define SPLINTER_BUILDER_H

namespace SPLINTER
{

/**
 * Pure abstract builder base class.
 * Used for building Functions
 */
template <class T>
class Builder
{
public:
    Builder() : _lambda(0.0) {}
    virtual ~Builder() {};

    virtual T build() const = 0;

    Builder<T> &lambda(double lambda)
    {
        if (lambda < 0)
        {
            throw Exception("Builder::lambda: Lambda must be non-negative.");
        }

        _lambda = lambda;
        return *this;
    }

private:
    double _lambda;
};

} // namespace SPLINTER

#endif // SPLINTER_BUILDER_H

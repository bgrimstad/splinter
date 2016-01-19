/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_RBFBUILDER_H
#define SPLINTER_RBFBUILDER_H

#include "builderbase.h"
#include "rbfnetwork.h"

namespace SPLINTER
{

class SPLINTER_API RBFNetwork::Builder : public BuilderBase<RBFNetwork>
{
public:
    Builder(const DataTable &data)
        :
        BuilderBase(data),
        _type(RBFType::THIN_PLATE_SPLINE),
        _normalized(false),
        _precondition(false)
    {}

    // Set build options
    Builder& type(RBFType type)
    {
        _type = type;
        return *this;
    }

    Builder& normalized(bool normalized)
    {
        _normalized = normalized;
        return *this;
    }

    Builder& precondition(bool precondition)
    {
        _precondition = precondition;
        return *this;
    }

    // Build RBFNetwork
    RBFNetwork build() const override;

private:
    Builder();

    DenseMatrix computePreconditionMatrix() const;

    RBFType _type;
    bool _normalized;
    bool _precondition;
};

} // namespace SPLINTER

#endif // SPLINTER_RBFBUILDER_H

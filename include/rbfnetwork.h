/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_RBFNETWORK_H
#define SPLINTER_RBFNETWORK_H

#include "definitions.h"
#include "datatable.h"
#include "linearfunction.h"
#include "rbf.h"
#include "memory"

namespace SPLINTER
{

/*
 * Class for radial basis function splines.
 * The RBF splines support scattered sampling, but their construction require
 * the solution of a fairly ill-conditioned linear system. This drawback may be
 * alleviated by applying a pre-conditioner to the linear system.
 */
class SPLINTER_API RBFNetwork : public LinearFunction<DenseVector, DenseMatrix>
{
public:
    RBFNetwork(const char *filename);
    RBFNetwork(const std::string &filename);
    RBFNetwork(const DataTable &samples, RBFType type);
    RBFNetwork(const DataTable &samples, RBFType type, bool normalized);

    virtual RBFNetwork* clone() const { return new RBFNetwork(*this); }

    // Currently overriding since evalBasisJacobian is not implemented
    DenseMatrix evalJacobian(DenseVector x) const override;

    // Should return dense types
    DenseVector evalBasis(DenseVector x) const override;
    DenseMatrix evalBasisJacobian(DenseVector x) const override {
        // TODO: implement
        return DenseMatrix::Zero(getNumCoefficients(), numVariables);
    };

    void save(const std::string &fileName) const override;

    std::string getDescription() const override;

private:
    DataTable samples;
    bool normalized, precondition;

    // Store the type so we can reconstruct the object when deserializing
    RBFType type;
    std::shared_ptr<RBF> fn;

    DenseMatrix computePreconditionMatrix() const;

    void load(const std::string &fileName) override;

    friend class Serializer;
    friend bool operator==(const RBFNetwork &lhs, const RBFNetwork &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_RBFNETWORK_H

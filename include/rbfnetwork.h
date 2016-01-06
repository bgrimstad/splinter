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
#include "function.h"
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
class SPLINTER_API RBFNetwork : public Function
{
public:
    RBFNetwork(const char *filename);
    RBFNetwork(const std::string filename);
    RBFNetwork(const DataTable &samples, RBFType type);
    RBFNetwork(const DataTable &samples, RBFType type, bool normalized);

    virtual RBFNetwork * clone() const { return new RBFNetwork(*this); }

    double eval(DenseVector x) const;
    double eval(std::vector<double> x) const;
    DenseMatrix evalJacobian(DenseVector x) const;
    DenseMatrix evalHessian(DenseVector x) const { DenseMatrix h(numVariables, numVariables); h.fill(0.0); return h; }; // TODO: implement
    //    std::vector<double> getDomainUpperBound() const;
    //    std::vector<double> getDomainLowerBound() const;

    void save(const std::string fileName) const override;

    const std::string getDescription() const override;

private:
    RBFNetwork();

    DataTable samples;
    bool normalized, precondition;
    unsigned int numSamples;

    // Store the type so we can reconstruct the object when deserializing
    RBFType type;
    std::shared_ptr<RBF> fn;

    DenseMatrix weights;

    DenseMatrix computePreconditionMatrix() const;

    double dist(std::vector<double> x, std::vector<double> y) const;
    double dist(DataPoint x, DataPoint y) const;
    bool dist_sort(DataPoint x, DataPoint y) const;

    void load(const std::string fileName) override;

    friend class Serializer;
    friend bool operator==(const RBFNetwork &lhs, const RBFNetwork &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_RBFNETWORK_H

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_BASIS_H
#define SPLINTER_BSPLINE_BASIS_H

#include "definitions.h"
#include "bspline_basis_1d.h"


namespace SPLINTER
{

class BSplineBasis
{
public:
    // Constructor
    BSplineBasis(std::vector<unsigned int> degrees, const std::vector<std::vector<double>> &knot_vectors);

    // Evaluation
    SparseVector eval(const DenseVector &x) const;
    DenseMatrix eval_basis_jacobian_old(const DenseVector &x) const; // Deprecated
    SparseMatrix eval_basis_jacobian(const DenseVector &x) const;
    SparseMatrix eval_basis_jacobian2(const DenseVector &x) const; // A bit slower than evaBasisJacobianOld()
    SparseMatrix eval_basis_hessian(const DenseVector &x) const;

    // Knot vector manipulation
    SparseMatrix refine_knots();
    SparseMatrix refine_knots_locally(const DenseVector &x);
    SparseMatrix decompose_to_bezier();
    SparseMatrix insert_knots(double tau, unsigned int dim, unsigned int multiplicity = 1);

    // Getters
    BSplineBasis1D get_single_basis(unsigned int dim);
    std::vector<std::vector<double>> get_knot_vectors() const;
    std::vector<double> get_knot_vector(int dim) const;

    std::vector<unsigned int> get_basis_degrees() const;
    unsigned int get_basis_degree(unsigned int dim) const;
    unsigned int get_num_basis_functions() const;
    unsigned int get_num_basis_functions(unsigned int dim) const;
    std::vector<unsigned int> get_num_basis_functions_target() const;

    unsigned int get_knot_multiplicity(unsigned int dim, double tau) const;

    /*
     * Returns the maximum number of supported basis functions at any point in the B-spline domain
     */
    unsigned int num_supported() const;

    bool inside_support(const DenseVector &x) const;
    std::vector<double> get_support_lower_bound() const;
    std::vector<double> get_support_upper_bound() const;

    // Support related
    SparseMatrix reduce_support(const std::vector<double> &lb, const std::vector<double> &ub);

private:
    std::vector<BSplineBasis1D> bases;
    unsigned int num_variables;

    friend bool operator==(const BSplineBasis &lhs, const BSplineBasis &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BASIS_H

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_BASIS_1D_H
#define SPLINTER_BSPLINE_BASIS_1D_H

#include "definitions.h"
#include "knot_vector.h"

namespace SPLINTER
{

class BSplineBasis1D
{
public:
    BSplineBasis1D(unsigned int degree, const std::vector<double> &knots);

    /**
     * Evaluation of basis functions
     */
    SparseVector eval(double x) const;
    SparseVector eval_derivative(double x, unsigned int r) const;
    SparseVector eval_first_derivative(double x) const; // TODO: Deprecated

    /**
     * Knot vector related
     */
    SparseMatrix refine_knots();
    SparseMatrix refine_knots_locally(double x);
    SparseMatrix decompose_to_bezier();
    SparseMatrix insert_knots(double tau, unsigned int multiplicity = 1);
    // bool insert_knots(SparseMatrix &A, std::vector<tuple<double,int>> newKnots); // Add knots at several locations

    unsigned int knot_multiplicity(double tau) const {
        // Return the number of repetitions of tau in the knot vector
        return knots.multiplicity(tau);
    }

    /**
     * Support related
     */
    double support_hack(double x) const;
    bool is_supported(double x) const {
        return knots.is_supported(x);
    }
    SparseMatrix reduce_support(double lb, double ub);

    /**
     * Getters
     */
    std::vector<double> get_knot_vector() const { return knots.get_values(); }
    unsigned int get_basis_degree() const { return degree; }
    unsigned int get_num_basis_functions() const;
    unsigned int get_num_basis_functions_target() const;

    /**
     * Index getters
     */
    std::vector<unsigned int> index_supported_basis_functions(double x) const;
    unsigned int index_longest_interval(const std::vector<double> &vec) const;

    /**
     * Setters
     */
    void set_num_basis_functions_target(unsigned int target)
    {
        target_num_basis_functions = std::max(degree+1, target);
    }

private:
    unsigned int degree;
    KnotVector knots;
    unsigned int target_num_basis_functions;

    // DeBoorCox algorithm for evaluating basis functions
    double de_boor_cox(double x, unsigned int i, unsigned int k) const;
    double de_boor_cox_coeff(double x, double x_min, double x_max) const;

    // Builds basis matrix for alternative evaluation of basis functions
    SparseMatrix build_basis_matrix(double x, unsigned int u, unsigned int k, bool diff = false) const;

    /**
     * Builds knot insertion matrix
     * Implements Oslo Algorithm 1 from Lyche and Moerken (2011). Spline methods draft.
     */
    SparseMatrix build_knot_insertion_matrix(const std::vector<double> &refined_knots) const;

    // Operators
    friend bool operator==(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs);
    friend bool operator!=(const BSplineBasis1D &lhs, const BSplineBasis1D &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_BASIS_1D_H

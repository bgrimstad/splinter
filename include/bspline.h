/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_BSPLINE_H
#define SPLINTER_BSPLINE_H

#include "function.h"
#include "bspline_basis.h"
#include "data_table.h"
#include "json_parser.h"


namespace SPLINTER
{

/**
 * Class that implements the multivariate tensor product B-spline
 */
class SPLINTER_API BSpline : public Function
{
public:

    // B-spline smoothing
    enum class Smoothing {
        NONE,       // No smoothing
        IDENTITY,   // Regularization term alpha*c'*I*c is added to OLS objective
        PSPLINE     // Smoothing term alpha*Delta(c,2) is added to OLS objective
    };

    /**
     * Builder class for construction by regression
     * Implemented in bspline_builder.*
     */
    class Builder;
//    enum class Smoothing;

    /**
     * Construct B-spline from knot vectors, control points, and basis degrees
     */
    BSpline(unsigned int dim_x,
            unsigned int dim_y,
            const std::vector<std::vector<double>> &knot_vectors,
            const std::vector<unsigned int> &degrees);

    BSpline(const std::vector<std::vector<double>> &control_points,
            const std::vector<std::vector<double>> &knot_vectors,
            const std::vector<unsigned int> &degrees);

    virtual BSpline* clone() const { return new BSpline(*this); }

    /**
     * Evaluation of B-spline
     */

    // Avoid name hiding
    using Function::eval;
    using Function::eval_jacobian;

    /**
     * Returns the (dimY) function values at x
     */
    std::vector<double> eval(const std::vector<double> &x) const override;
    DenseVector eval(const DenseVector &x) const override;

    /**
     * Returns the (dimY x dimX) Jacobian evaluated at x
     */
    DenseMatrix eval_jacobian(const DenseVector &x) const override;

    /**
     * Returns the (dimY x dimX x dimX) Hessian evaluated at x
     */
    std::vector<std::vector<std::vector<double>>> eval_hessian(const std::vector<double> &x) const;

    // Evaluation of B-spline basis functions
    SparseVector eval_basis(const DenseVector &x) const;
    SparseMatrix eval_basis_jacobian(const DenseVector &x) const;

    /**
     * Getters
     */
    DenseMatrix get_control_points() const
    {
        return control_points;
    }

    unsigned int get_num_control_points() const
    {
        return (unsigned int) control_points.rows();
    }

    std::vector<unsigned int> get_num_basis_functions_per_variable() const;

    unsigned int get_num_basis_functions() const
    {
        return basis.get_num_basis_functions();
    }

//    DenseMatrix get_control_points() const;
    DenseMatrix get_knot_averages() const {
        return compute_knot_averages();
    };
    std::vector< std::vector<double>> get_knot_vectors() const;
    std::vector<unsigned int> get_basis_degrees() const;
    std::vector<double> get_domain_upper_bound() const;
    std::vector<double> get_domain_lower_bound() const;

    /**
     * Setters
     */
    void set_control_points(const DenseMatrix &new_control_points);

    /**
     * Manipulations to B-spline
     */

    // Linear transformation of control points (B-spline has affine invariance)
    void linear_transform(const SparseMatrix &A);

    // Fit B-spline to sample data
    BSpline& fit(const DataTable &data, Smoothing smoothing = Smoothing::NONE, double alpha = .1,
                 std::vector<double> weights = std::vector<double>());

    // Reduce support of B-spline
    void reduce_support(const std::vector<double> &lb, const std::vector<double> &ub,
                        bool regularize_knot_vectors = true);

    // Perform global knot refinement (all knots in one shabang)
    void global_knot_refinement();

    // Perform a local knot refinement at x
    void local_knot_refinement(const DenseVector &x);

    // Decompose B-spline to Bezier form
    void decompose_to_bezier();

    // Insert a knot until desired knot multiplicity is obtained
    void insert_knots(double tau, unsigned int dim, unsigned int multiplicity = 1);

    /**
     * Save and load
     */
    void to_json(const std::string &filename) const;

    static BSpline from_json(const std::string &filename);

    /**
     * Helper functions
     */
    void check_control_points() const;
    std::string get_description() const override;

protected:
    BSpline();

    BSplineBasis basis;

    /*
     * B-spline control points in R^(m x n),
     * where m = num_basis_functions and n = dim_y.
     * Each row is a control point.
     */
    DenseMatrix control_points;

    // Control point computations
    DenseMatrix compute_knot_averages() const;

private:
    // Domain reduction
    void regularize_knot_vectors(const std::vector<double> &lb, const std::vector<double> &ub);
    bool remove_unsupported_basis_functions(const std::vector<double> &lb, const std::vector<double> &ub);

    // Helper functions
    bool is_supported(const DenseVector &x) const;

    friend bool operator==(const BSpline &lhs, const BSpline &rhs);
};

} // namespace SPLINTER

#endif // SPLINTER_BSPLINE_H

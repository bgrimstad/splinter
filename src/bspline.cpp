/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline.h"
#include "bspline_basis.h"
#include "kronecker_product.h"
#include "unsupported/Eigen/KroneckerProduct"
#include <linear_solvers.h>
#include <iostream>
#include <utilities.h>
#include "bspline_utils.h"
#include "knot_builders.h"


namespace SPLINTER
{

/**
 * Constructors for multivariate B-spline using explicit data
 */
BSpline::BSpline(const std::vector<unsigned int> &degrees,
                 const std::vector<std::vector<double>> &knot_vectors,
                 unsigned int dim_y)
    : Function(degrees.size(), dim_y),
      basis(BSplineBasis(degrees, knot_vectors)),
      control_points(DenseMatrix::Zero(basis.get_num_basis_functions(), dim_y))
{
    check_control_points();
}

BSpline::BSpline(const std::vector<unsigned int> &degrees,
                 const std::vector<std::vector<double>> &knot_vectors,
                 const std::vector<std::vector<double>> &control_points)
    : Function(knot_vectors.size(), control_points.at(0).size()),
      basis(BSplineBasis(degrees, knot_vectors)),
      control_points(std_to_eig_mat(control_points))
{
    check_control_points();
}

/**
 * Returns the function value at x
 */
std::vector<double> BSpline::eval(const std::vector<double> &x) const {
    return eig_to_std_vec(eval(std_to_eig_vec(x)));
}

DenseVector BSpline::eval(const DenseVector &x) const
{
    check_input(x);
    DenseVector res = control_points.transpose()* eval_basis(x);
    return res;
}

/**
 * Returns the (dimY x dimX) Jacobian evaluated at x
 */
DenseMatrix BSpline::eval_jacobian(const DenseVector &x) const
{
    check_input(x);
    return control_points.transpose()* eval_basis_jacobian(x);
}

/**
 * Returns the Hessian evaluated at x. The Hessian is a (dim_y x dim_x x dim_x) tensor.
 */
std::vector<std::vector<std::vector<double>>> BSpline::eval_hessian(const std::vector<double> &x) const
{
    DenseVector eigX = std_to_eig_vec(x);
    check_input(eigX);

    std::vector<std::vector<std::vector<double>>> hessian;

    DenseMatrix identity = DenseMatrix::Identity(dim_x, dim_x);

    DenseMatrix cp_copy = DenseMatrix(control_points);

    for (size_t i = 0; i < dim_y; ++i)
    {
        DenseMatrix H = DenseMatrix::Zero(1, 1);
        DenseMatrix cp = cp_copy.col(i);
        DenseMatrix caug = kroneckerProduct(identity, cp.transpose());
        DenseMatrix DB = basis.eval_basis_hessian(eigX);

        H = caug*DB;

//        std::cout << cp << std::endl;

        // Fill in upper triangular of Hessian
        for (size_t j = 0; j < dim_x; ++j)
            for (size_t k = j+1; k < dim_x; ++k)
                H(j, k) = H(k, j);

        hessian.push_back(eig_to_std_mat(H));
    }

    return hessian;
}

/**
 * Evaluation of B-spline basis functions
 */
SparseVector BSpline::eval_basis(const DenseVector &x) const
{
#ifndef NDEBUG
    if (!is_supported(x))
        std::cout << "BSpline::eval_basis: Evaluation at point outside of support." << std::endl;
#endif // NDEBUG

    return basis.eval(x);
}

/**
 * Evaluation of B-spline basis Jacobian
 */
SparseMatrix BSpline::eval_basis_jacobian(const DenseVector &x) const
{
#ifndef NDEBUG
    if (!is_supported(x))
        std::cout << "BSpline::eval_basis_jacobian: Evaluation at point outside of support." << std::endl;
#endif // NDEBUG

    //SparseMatrix Bi = basis.eval_basis_jacobian(x);       // Sparse Jacobian implementation
    //SparseMatrix Bi = basis.eval_basis_jacobian2(x);  // Sparse Jacobian implementation
    DenseMatrix Bi = basis.eval_basis_jacobian_old(x);  // Old Jacobian implementation

    return Bi.sparseView();
}

std::vector<unsigned int> BSpline::get_num_basis_functions_per_variable() const
{
    std::vector<unsigned int> ret;
    for (unsigned int i = 0; i < dim_x; i++)
        ret.push_back(basis.get_num_basis_functions(i));
    return ret;
}

std::vector<std::vector<double>> BSpline::get_knot_vectors() const
{
    return basis.get_knot_vectors();
}

std::vector<unsigned int> BSpline::get_basis_degrees() const
{
    return basis.get_basis_degrees();
}

std::vector<double> BSpline::get_domain_upper_bound() const
{
    return basis.get_support_upper_bound();
}

std::vector<double> BSpline::get_domain_lower_bound() const
{
    return basis.get_support_lower_bound();
}

void BSpline::set_control_points(const DenseMatrix &new_control_points)
{
    if (new_control_points.rows() != get_num_basis_functions())
        throw Exception("BSpline::set_control_points: Incompatible size of coefficient vector. " +
                                std::to_string(new_control_points.rows()) + " not equal to " +
                                std::to_string(get_num_basis_functions()) + "!");

    this->control_points = new_control_points;
    check_control_points();
}

void BSpline::linear_transform(const SparseMatrix &A)
{
    if (A.cols() != control_points.rows())
        throw Exception("BSpline::linear_transform: Incompatible size of linear transformation matrix.");
    control_points = A*control_points;
}

BSpline& BSpline::fit(const DataTable &data, Smoothing smoothing, double alpha, std::vector<double> weights)
{
    if (data.get_dim_x() != get_dim_x()) {
        throw Exception("BSpline::fit: Expected " + std::to_string(get_dim_x()) + " input variables.");
    }

    if (data.get_dim_y() != get_dim_y()) {
        throw Exception("BSpline::fit: Expected " + std::to_string(get_dim_y()) + " output variables.");
    }

    if (alpha < 0) {
        throw Exception("BSpline::fit: alpha must be non-negative.");
    }

    if (weights.size() > 0 && data.get_num_samples() != weights.size()) {
        throw Exception("BSpline::fit: number of weights must equal number of data points.");
    }

#ifndef NDEBUG
    if (!data.is_grid_complete())
        std::cout << "BSpline::fit: Fitting B-spline to scattered data points (irregular grid of points)." << std::endl;
#endif // NDEBUG

    // Compute control points from samples and update B-spline
    auto coefficients = compute_control_points(*this, data, smoothing, alpha, weights);
    set_control_points(coefficients);

    return *this;
}

void BSpline::check_control_points() const
{
    if (control_points.cols() != get_dim_y())
        throw Exception("BSpline::check_control_points: Inconsistent number of columns of control points matrix.");
    if (control_points.rows() != get_num_basis_functions())
        throw Exception("BSpline::check_control_points: Inconsistent number of rows of control points matrix.");
}

bool BSpline::is_supported(const DenseVector &x) const
{
    return basis.inside_support(x);
}

void BSpline::reduce_support(const std::vector<double> &lb, const std::vector<double> &ub, bool regularize_knot_vectors)
{
    if (lb.size() != dim_x || ub.size() != dim_x)
        throw Exception("BSpline::reduce_support: Inconsistent vector sizes!");

    std::vector<double> sl = basis.get_support_lower_bound();
    std::vector<double> su = basis.get_support_upper_bound();

    for (unsigned int dim = 0; dim < dim_x; dim++)
    {
        // Check if new domain is empty
        if (ub.at(dim) <= lb.at(dim) || lb.at(dim) >= su.at(dim) || ub.at(dim) <= sl.at(dim))
            throw Exception("BSpline::reduce_support: Cannot reduce B-spline domain to empty set!");

        // Check if new domain is a strict subset
        if (su.at(dim) < ub.at(dim) || sl.at(dim) > lb.at(dim))
            throw Exception("BSpline::reduce_support: Cannot expand B-spline domain!");

        // Tightest possible
        sl.at(dim) = lb.at(dim);
        su.at(dim) = ub.at(dim);
    }

    if (regularize_knot_vectors)
    {
        this->regularize_knot_vectors(sl, su);
    }

    // Remove knots and control points that are unsupported with the new bounds
    if (!remove_unsupported_basis_functions(sl, su))
    {
        throw Exception("BSpline::reduce_support: Failed to remove unsupported basis functions!");
    }
}

void BSpline::global_knot_refinement()
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.refine_knots();

    // Update control points
    linear_transform(A);
}

void BSpline::local_knot_refinement(const DenseVector &x)
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.refine_knots_locally(x);

    // Update control points
    linear_transform(A);
}

void BSpline::decompose_to_bezier()
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.decompose_to_bezier();

    // Update control points
    linear_transform(A);
}

// Computes knot averages: assumes that basis is initialized!
DenseMatrix BSpline::compute_knot_averages() const
{
    // Calculate knot averages for each knot vector
    std::vector<DenseVector> mu_vectors;
    for (unsigned int i = 0; i < dim_x; i++)
    {
        std::vector<double> knots = basis.get_knot_vector(i);
        DenseVector mu = DenseVector::Zero(basis.get_num_basis_functions(i));

        for (unsigned int j = 0; j < basis.get_num_basis_functions(i); j++)
        {
            double knotAvg = 0;
            for (unsigned int k = j+1; k <= j+ basis.get_basis_degree(i); k++)
            {
                knotAvg += knots.at(k);
            }
            mu(j) = knotAvg/ basis.get_basis_degree(i);
        }
        mu_vectors.push_back(mu);
    }

    // Calculate vectors of ones (with same length as corresponding knot average vector)
    std::vector<DenseVector> knot_ones;
    for (unsigned int i = 0; i < dim_x; i++)
        knot_ones.push_back(DenseVector::Ones(mu_vectors.at(i).rows()));

    // Fill knot average matrix one column at the time
    DenseMatrix knot_averages = DenseMatrix::Zero(basis.get_num_basis_functions(), dim_x);

    for (unsigned int i = 0; i < dim_x; i++)
    {
        DenseMatrix mu_ext(1,1); mu_ext(0,0) = 1;
        for (unsigned int j = 0; j < dim_x; j++)
        {
            DenseMatrix temp = mu_ext;
            if (i == j)
                mu_ext = Eigen::kroneckerProduct(temp, mu_vectors.at(j));
            else
                mu_ext = Eigen::kroneckerProduct(temp, knot_ones.at(j));
        }
        if (mu_ext.rows() != basis.get_num_basis_functions())
            throw Exception("BSpline::compute_knot_averages: Incompatible size of knot average matrix.");
        knot_averages.block(0, i, basis.get_num_basis_functions(), 1) = mu_ext;
    }

    return knot_averages;
}

void BSpline::insert_knots(double tau, unsigned int dim, unsigned int multiplicity)
{
    // Insert knots and compute knot insertion matrix
    SparseMatrix A = basis.insert_knots(tau, dim, multiplicity);

    // Update control points
    linear_transform(A);
}

void BSpline::regularize_knot_vectors(const std::vector<double> &lb, const std::vector<double> &ub)
{
    // Add and remove controlpoints and knots to make the B-spline p-regular with support [lb, ub]
    if (!(lb.size() == dim_x && ub.size() == dim_x))
        throw Exception("BSpline::regularize_knot_vectors: Inconsistent vector sizes.");

    for (unsigned int dim = 0; dim < dim_x; dim++)
    {
        unsigned int multiplicity_target = basis.get_basis_degree(dim) + 1;

        // Inserting many knots at the time (to save number of B-spline coefficient calculations)
        // NOTE: This method generates knot insertion matrices with more nonzero elements than
        // the method that inserts one knot at the time. This causes the preallocation of
        // kronecker product matrices to become too small and the speed deteriorates drastically
        // in higher dimensions because reallocation is necessary. This can be prevented by
        // computing the number of nonzeros when preallocating memory (see my_kronecker_product).
        int num_knots_lb = multiplicity_target - basis.get_knot_multiplicity(dim, lb.at(dim));
        if (num_knots_lb > 0)
        {
            insert_knots(lb.at(dim), dim, num_knots_lb);
        }

        int num_knots_ub = multiplicity_target - basis.get_knot_multiplicity(dim, ub.at(dim));
        if (num_knots_ub > 0)
        {
            insert_knots(ub.at(dim), dim, num_knots_ub);
        }
    }
}

bool BSpline::remove_unsupported_basis_functions(const std::vector<double> &lb, const std::vector<double> &ub)
{
    if (lb.size() != dim_x || ub.size() != dim_x)
        throw Exception("BSpline::remove_unsupported_basis_functions: Incompatible dimension of domain bounds.");

    SparseMatrix A = basis.reduce_support(lb, ub);

    // Remove unsupported control points (basis functions)
    linear_transform(A);

    return true;
}

void BSpline::to_json(const std::string &filename) const {
    SPLINTER::bspline_to_json(*this, filename);
}

BSpline BSpline::from_json(const std::string &filename) {
    return SPLINTER::bspline_from_json(filename);
}

std::string BSpline::get_description() const
{
    std::string description("BSpline of degree");
    auto degrees = get_basis_degrees();
    // See if all degrees are the same.
    bool equal = true;
    for (size_t i = 1; i < degrees.size(); ++i)
    {
        equal = equal && (degrees.at(i) == degrees.at(i-1));
    }

    if(equal)
    {
        description.append(" ");
        description.append(std::to_string(degrees.at(0)));
    }
    else
    {
        description.append("s (");
        for (size_t i = 0; i < degrees.size(); ++i)
        {
            description.append(std::to_string(degrees.at(i)));
            if (i + 1 < degrees.size())
            {
                description.append(", ");
            }
        }
        description.append(")");
    }

    return description;
}

} // namespace SPLINTER

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
#include <serializer.h>
#include <iostream>
#include <utilities.h>


namespace SPLINTER
{

/**
 * Constructor used by serializer
 */
BSpline::BSpline()
    : Function()
{}

/**
 * Constructors for multivariate B-spline using explicit data
 */
BSpline::BSpline(
        unsigned int dimX,
        unsigned int dimY,
        const std::vector<std::vector<double>> &knotVectors,
        const std::vector<unsigned int> &degrees)
    : Function(dimX, dimY),
      basis(BSplineBasis(knotVectors, degrees)),
      control_points(DenseMatrix::Zero(basis.getNumBasisFunctions(), dimY))
{
    check_control_points();
}

BSpline::BSpline(const std::vector<std::vector<double>> &controlPoints,
                 const std::vector<std::vector<double>> &knotVectors,
                 const std::vector<unsigned int> &degrees)
    : Function(knotVectors.size(), controlPoints.at(0).size()),
      basis(BSplineBasis(knotVectors, degrees)),
      control_points(stdVecVecToEigMat(controlPoints))
{
    check_control_points();
}

/**
 * Construct B-spline from saved data
 */
BSpline::BSpline(const char *fileName)
    : BSpline(std::string(fileName))
{
}

BSpline::BSpline(const std::string &fileName)
    : Function()
{
    load(fileName);
}

/**
 * Returns the function value at x
 */
std::vector<double> BSpline::eval(const std::vector<double> &x) const {
    return eigToStdVec(eval(stdToEigVec(x)));
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

/*
 * Returns the Hessian evaluated at x. The Hessian is a (dim_y x dim_x x dim_x) tensor.
 */
std::vector<std::vector<std::vector<double>>> BSpline::eval_hessian(const std::vector<double> &x) const
{
    DenseVector eigX = stdToEigVec(x);
    check_input(eigX);

    std::vector<std::vector<std::vector<double>>> hessian;

    DenseMatrix identity = DenseMatrix::Identity(dim_x, dim_x);

    DenseMatrix cpCopy = DenseMatrix(control_points);

    for (size_t i = 0; i < dim_y; ++i)
    {
        DenseMatrix H = DenseMatrix::Zero(1, 1);
        DenseMatrix cp = cpCopy.col(i);
        DenseMatrix caug = kroneckerProduct(identity, cp.transpose());
        DenseMatrix DB = basis.evalBasisHessian(eigX);

        H = caug*DB;

//        std::cout << cp << std::endl;

        // Fill in upper triangular of Hessian
        for (size_t j = 0; j < dim_x; ++j)
            for (size_t k = j+1; k < dim_x; ++k)
                H(j, k) = H(k, j);

        hessian.push_back(eigMatToStdVecVec(H));
    }

    return hessian;
}

// Evaluation of B-spline basis functions
SparseVector BSpline::eval_basis(const DenseVector &x) const
{
#ifndef NDEBUG
    if (!isSupported(x))
        std::cout << "BSpline::evalBasis: Evaluation at point outside of support." << std::endl;
#endif // NDEBUG

    return basis.eval(x);
}

SparseMatrix BSpline::eval_basis_jacobian(const DenseVector &x) const
{
#ifndef NDEBUG
    if (!isSupported(x))
        std::cout << "BSpline::eval_basis_jacobian: Evaluation at point outside of support." << std::endl;
#endif // NDEBUG

    //SparseMatrix Bi = basis.eval_basis_jacobian(x);       // Sparse Jacobian implementation
    //SparseMatrix Bi = basis.evalBasisJacobian2(x);  // Sparse Jacobian implementation
    DenseMatrix Bi = basis.evalBasisJacobianOld(x);  // Old Jacobian implementation

    return Bi.sparseView();
}

std::vector<unsigned int> BSpline::get_num_basis_functions_per_variable() const
{
    std::vector<unsigned int> ret;
    for (unsigned int i = 0; i < dim_x; i++)
        ret.push_back(basis.getNumBasisFunctions(i));
    return ret;
}

std::vector<std::vector<double>> BSpline::get_knot_vectors() const
{
    return basis.getKnotVectors();
}

std::vector<unsigned int> BSpline::get_basis_degrees() const
{
    return basis.getBasisDegrees();
}

std::vector<double> BSpline::get_domain_upper_bound() const
{
    return basis.getSupportUpperBound();
}

std::vector<double> BSpline::get_domain_lower_bound() const
{
    return basis.getSupportLowerBound();
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

void BSpline::check_control_points() const
{
    if (control_points.cols() != get_dim_y())
        throw Exception("BSpline::check_control_points: Inconsistent number of columns of control points matrix.");
    if (control_points.rows() != get_num_basis_functions())
        throw Exception("BSpline::check_control_points: Inconsistent number of rows of control points matrix.");
}

bool BSpline::is_supported(const DenseVector &x) const
{
    return basis.insideSupport(x);
}

void BSpline::reduce_support(const std::vector<double> &lb, const std::vector<double> &ub, bool regularize_knot_vectors)
{
    if (lb.size() != dim_x || ub.size() != dim_x)
        throw Exception("BSpline::reduce_support: Inconsistent vector sizes!");

    std::vector<double> sl = basis.getSupportLowerBound();
    std::vector<double> su = basis.getSupportUpperBound();

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
    SparseMatrix A = basis.refineKnots();

    // Update control points
    linear_transform(A);
}

void BSpline::local_knot_refinement(const DenseVector &x)
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.refineKnotsLocally(x);

    // Update control points
    linear_transform(A);
}

void BSpline::decompose_to_bezier()
{
    // Compute knot insertion matrix
    SparseMatrix A = basis.decomposeToBezierForm();

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
        std::vector<double> knots = basis.getKnotVector(i);
        DenseVector mu = DenseVector::Zero(basis.getNumBasisFunctions(i));

        for (unsigned int j = 0; j < basis.getNumBasisFunctions(i); j++)
        {
            double knotAvg = 0;
            for (unsigned int k = j+1; k <= j+basis.getBasisDegree(i); k++)
            {
                knotAvg += knots.at(k);
            }
            mu(j) = knotAvg/basis.getBasisDegree(i);
        }
        mu_vectors.push_back(mu);
    }

    // Calculate vectors of ones (with same length as corresponding knot average vector)
    std::vector<DenseVector> knotOnes;
    for (unsigned int i = 0; i < dim_x; i++)
        knotOnes.push_back(DenseVector::Ones(mu_vectors.at(i).rows()));

    // Fill knot average matrix one column at the time
    DenseMatrix knot_averages = DenseMatrix::Zero(basis.getNumBasisFunctions(), dim_x);

    for (unsigned int i = 0; i < dim_x; i++)
    {
        DenseMatrix mu_ext(1,1); mu_ext(0,0) = 1;
        for (unsigned int j = 0; j < dim_x; j++)
        {
            DenseMatrix temp = mu_ext;
            if (i == j)
                mu_ext = Eigen::kroneckerProduct(temp, mu_vectors.at(j));
            else
                mu_ext = Eigen::kroneckerProduct(temp, knotOnes.at(j));
        }
        if (mu_ext.rows() != basis.getNumBasisFunctions())
            throw Exception("BSpline::compute_knot_averages: Incompatible size of knot average matrix.");
        knot_averages.block(0, i, basis.getNumBasisFunctions(), 1) = mu_ext;
    }

    return knot_averages;
}

void BSpline::insert_knots(double tau, unsigned int dim, unsigned int multiplicity)
{
    // Insert knots and compute knot insertion matrix
    SparseMatrix A = basis.insertKnots(tau, dim, multiplicity);

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
        unsigned int multiplicityTarget = basis.getBasisDegree(dim) + 1;

        // Inserting many knots at the time (to save number of B-spline coefficient calculations)
        // NOTE: This method generates knot insertion matrices with more nonzero elements than
        // the method that inserts one knot at the time. This causes the preallocation of
        // kronecker product matrices to become too small and the speed deteriorates drastically
        // in higher dimensions because reallocation is necessary. This can be prevented by
        // computing the number of nonzeros when preallocating memory (see myKroneckerProduct).
        int numKnotsLB = multiplicityTarget - basis.getKnotMultiplicity(dim, lb.at(dim));
        if (numKnotsLB > 0)
        {
            insert_knots(lb.at(dim), dim, numKnotsLB);
        }

        int numKnotsUB = multiplicityTarget - basis.getKnotMultiplicity(dim, ub.at(dim));
        if (numKnotsUB > 0)
        {
            insert_knots(ub.at(dim), dim, numKnotsUB);
        }
    }
}

bool BSpline::remove_unsupported_basis_functions(const std::vector<double> &lb, const std::vector<double> &ub)
{
    if (lb.size() != dim_x || ub.size() != dim_x)
        throw Exception("BSpline::remove_unsupported_basis_functions: Incompatible dimension of domain bounds.");

    SparseMatrix A = basis.reduceSupport(lb, ub);

    // Remove unsupported control points (basis functions)
    linear_transform(A);

    return true;
}

void BSpline::save(const std::string &fileName) const
{
    Serializer s;
    s.serialize(*this);
    s.saveToFile(fileName);
}

void BSpline::load(const std::string &fileName)
{
    Serializer s(fileName);
    s.deserialize(*this);
}

void BSpline::save_to_json(const std::string &filename) const {
    SPLINTER::save_to_json(*this, filename);
}

BSpline BSpline::load_from_json(const std::string &filename) {
    return SPLINTER::load_from_json(filename);
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

/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <bspline_basis_1d.h>
#include <knot_builders.h>
#include <algorithm>
#include <utilities.h>
#include <iostream>

namespace SPLINTER
{

BSplineBasis1D::BSplineBasis1D(unsigned int degree, const std::vector<double> &knots)
    : degree(degree),
      knots(KnotVector(knots)),
      target_num_basis_functions((degree+1)+2*degree+1) // Minimum p+1
{
    // Check that knot vector is (p+1)-regular
    if (!this->knots.is_regular(degree))
        throw Exception("BSplineBasis1D::BSplineBasis1D: Knot vector is not regular.");
}

SparseVector BSplineBasis1D::eval(double x) const
{
    SparseVector values(get_num_basis_functions());

    if (!is_supported(x))
        return values;

    x = support_hack(x);

    auto index_supported = index_supported_basis_functions(x);

    values.reserve(index_supported.size());

    // Evaluate nonzero basis functions
    for (auto i : index_supported)
    {
        double val = de_boor_cox(x, i, degree);
        if (fabs(val) > 1e-12)
            values.insert(i) = val;
    }

    // Alternative evaluation using basis matrix
//    int knotIndex = indexHalfopenInterval(x); // knot index

//    SparseMatrix basisvalues2 = buildBsplineMatrix(x, knotIndex, 1);
//    for (int i = 2; i <= basisDegree; i++)
//    {
//        SparseMatrix Ri = buildBsplineMatrix(x, knotIndex, i);
//        basisvalues2 = basisvalues2*Ri;
//    }
//    basisvalues2.makeCompressed();

//    assert(basisvalues2.rows() == 1);
//    assert(basisvalues2.cols() == basisDegree + 1);

    return values;
}

SparseVector BSplineBasis1D::eval_derivative(double x, unsigned int r) const
{
    // Evaluate rth derivative of basis functions at x
    // Returns vector [D^(r)B_(u-p,p)(x) ... D^(r)B_(u,p)(x)]
    // where u is the knot index and p is the degree
    auto p = degree;

    // Continuity requirement
    if (p < r)
    {
        // Return zero-gradient
        SparseVector DB(get_num_basis_functions());
        return DB;
    }

    // TODO: Check for knot multiplicity here!

    x = support_hack(x);

    unsigned int knot_index = knots.index_interval(x);

    // Algorithm 3.18 from Lyche and Moerken (2011)
    SparseMatrix B(1,1);
    B.insert(0,0) = 1;

    for (unsigned int i = 1; i <= p-r; i++)
    {
        SparseMatrix R = build_basis_matrix(x, knot_index, i);
        B = B*R;
    }

    for (unsigned int i = p-r+1; i <= p; i++)
    {
        SparseMatrix DR = build_basis_matrix(x, knot_index, i, true);
        B = B*DR;
    }
    double factorial = std::tgamma(p+1)/std::tgamma(p-r+1);
    B = B*factorial;

    if (B.cols() != p+1)
        throw Exception("BSplineBasis1D::eval_derivative: Wrong number of columns of B matrix.");

    // From row vector to extended column vector
    SparseVector DB(get_num_basis_functions());
    DB.reserve(p+1);
    unsigned int i = knot_index-p; // First insertion index

    if (i < 0)
        throw Exception("BSplineBasis1D::eval_derivative: negative insertion index!");

    for (int k = 0; k < B.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(B, k); it; ++it)
        {
            DB.insert(i+it.col()) = it.value();
        }
    }

    return DB;
}

// Old implementation of first derivative of basis functions
SparseVector BSplineBasis1D::eval_first_derivative(double x) const
{
    SparseVector values(get_num_basis_functions());

    x = support_hack(x);

    auto supportedBasisFunctions = index_supported_basis_functions(x);

    for (auto i : supportedBasisFunctions)
    {
        // Differentiate basis function
        // Equation 3.35 in Lyche & Moerken (2011)
        double b1 = de_boor_cox(x, i, degree - 1);
        double b2 = de_boor_cox(x, i + 1, degree - 1);

        double t11 = knots.at(i);
        double t12 = knots.at(i+degree);
        double t21 = knots.at(i+1);
        double t22 = knots.at(i+degree+1);

        (t12 == t11) ? b1 = 0 : b1 = b1/(t12-t11);
        (t22 == t21) ? b2 = 0 : b2 = b2/(t22-t21);

        values.insert(i) = degree*(b1 - b2);
    }

    return values;
}

/* 
 * Used to evaluate basis functions - alternative to the recursive deBoorCox
 * Builds B-spline matrix R_k in R^(k,k+1), or, if diff = true, the differentiated basis matrix DR_k in R^(k,k+1).
 */ 
SparseMatrix BSplineBasis1D::build_basis_matrix(double x, unsigned int u, unsigned int k, bool diff) const
{
    if (!(k >= 1 && k <= get_basis_degree())) {
        throw Exception("BSplineBasis1D::build_basis_matrix: Incorrect input parameters!");
    }

//    assert(u >= basisDegree + 1);
//    assert(u < ks.size() - basisDegree);

    auto rows = k;
    auto cols = k+1;
    SparseMatrix R(rows, cols);
    R.reserve(Eigen::VectorXi::Constant(cols, 2));

    for (unsigned int i = 0; i < rows; i++)
    {
        double dk = knots.at(u+1+i) - knots.at(u+1+i-k);
        if (dk == 0)
        {
            continue;
        }
        else
        {
            if (diff)
            {
                // Insert diagonal element
                R.insert(i,i) = -1/dk;

                // Insert super-diagonal element
                R.insert(i,i+1) = 1/dk;
            }
            else
            {
                // Insert diagonal element
                double a = (knots.at(u+1+i) - x)/dk;
                if (!assert_near(a, .0))
                    R.insert(i,i) = a;

                // Insert super-diagonal element
                double b = (x - knots.at(u+1+i-k))/dk;
                if (!assert_near(b, .0))
                    R.insert(i,i+1) = b;
            }
        }
    }

    R.makeCompressed();

    return R;
}

double BSplineBasis1D::de_boor_cox(double x, unsigned int i, unsigned int k) const
{
    if (k == 0)
    {
        if ((knots.at(i) <= x) && (x < knots.at(i+1)))
            return 1;
        else
            return 0;
    }
    else
    {
        double s1,s2,r1,r2;

        s1 = de_boor_cox_coeff(x, knots.at(i), knots.at(i + k));
        s2 = de_boor_cox_coeff(x, knots.at(i + 1), knots.at(i + k + 1));

        r1 = de_boor_cox(x, i, k - 1);
        r2 = de_boor_cox(x, i + 1, k - 1);

        return s1*r1 + (1-s2)*r2;
    }
}

double BSplineBasis1D::de_boor_cox_coeff(double x, double x_min, double x_max) const
{
    if (x_min < x_max && x_min <= x && x <= x_max)
        return (x - x_min)/(x_max - x_min);
    return 0;
}

// Insert knots and compute knot insertion matrix (to update control points)
SparseMatrix BSplineBasis1D::insert_knots(double tau, unsigned int multiplicity)
{
    if (!is_supported(tau))
        throw Exception("BSplineBasis1D::insert_knots: Cannot insert knot outside domain!");

    if (knot_multiplicity(tau) + multiplicity > degree + 1)
        throw Exception("BSplineBasis1D::insert_knots: Knot multiplicity is too high!");

    // New knot vector
    int index = knots.index_interval(tau);

    std::vector<double> extended_knots = knots.get_values();
    for (unsigned int i = 0; i < multiplicity; i++)
        extended_knots.insert(extended_knots.begin()+index+1, tau);

    if (!KnotVector(extended_knots).is_regular(degree))
        throw Exception("BSplineBasis1D::insert_knots: New knot vector is not regular!");

    // Return knot insertion matrix
    SparseMatrix A = build_knot_insertion_matrix(extended_knots);

    // Update knots
    knots = KnotVector(extended_knots);

    return A;
}

SparseMatrix BSplineBasis1D::refine_knots()
{
    // Build refine knot vector
    std::vector<double> refinedKnots = knots.get_values();

    unsigned int targetNumKnots = target_num_basis_functions + degree + 1;
    while (refinedKnots.size() < targetNumKnots)
    {
        int index = index_longest_interval(refinedKnots);
        double newKnot = (refinedKnots.at(index) + refinedKnots.at(index+1))/2.0;
        refinedKnots.insert(std::lower_bound(refinedKnots.begin(), refinedKnots.end(), newKnot), newKnot);
    }

    if (!KnotVector(refinedKnots).is_regular(degree))
        throw Exception("BSplineBasis1D::refine_knots: New knot vector is not regular!");

    if (!knots.is_refinement(refinedKnots))
        throw Exception("BSplineBasis1D::refine_knots: New knot vector is not a proper refinement!");

    // Return knot insertion matrix
    SparseMatrix A = build_knot_insertion_matrix(refinedKnots);

    // Update knots
    knots = KnotVector(refinedKnots);

    return A;
}

SparseMatrix BSplineBasis1D::refine_knots_locally(double x)
{
    if (!is_supported(x))
        throw Exception("BSplineBasis1D::refine_knots_locally: Cannot refine outside support!");

    if (get_num_basis_functions() >= get_num_basis_functions_target()
            || assert_near(knots.front(), knots.back()))
    {
        unsigned int n = get_num_basis_functions();
        DenseMatrix A = DenseMatrix::Identity(n, n);
        return A.sparseView();
    }

    // Refined knot vector
    std::vector<double> refinedKnots = knots.get_values();

    auto upper = std::lower_bound(refinedKnots.begin(), refinedKnots.end(), x);

    // Check left boundary
    if (upper == refinedKnots.begin())
        std::advance(upper, degree+1);

    // Get previous iterator
    auto lower = std::prev(upper);

    // Do not insert if upper and lower bounding knot are close
    if (assert_near(*upper, *lower))
    {
        unsigned int n = get_num_basis_functions();
        DenseMatrix A = DenseMatrix::Identity(n,n);
        return A.sparseView();
    }

    // Insert knot at x
    double insertVal = x;

    // Adjust x if it is on or close to a knot
    if (knot_multiplicity(x) > 0
            || assert_near(*upper, x, 1e-6, 1e-6)
            || assert_near(*lower, x, 1e-6, 1e-6))
    {
        insertVal = (*upper + *lower)/2.0;
    }

    // Insert new knot
    refinedKnots.insert(upper, insertVal);

    if (!KnotVector(refinedKnots).is_regular(degree))
        throw Exception("BSplineBasis1D::refine_knots_locally: New knot vector is not regular!");

    if (!knots.is_refinement(refinedKnots))
        throw Exception("BSplineBasis1D::refine_knots_locally: New knot vector is not a proper refinement!");

    // Build knot insertion matrix
    SparseMatrix A = build_knot_insertion_matrix(refinedKnots);

    // Update knots
    knots = KnotVector(refinedKnots);

    return A;
}

SparseMatrix BSplineBasis1D::decompose_to_bezier()
{
    // Build refine knot vector
    std::vector<double> refined_knots = knots.get_values();

    // Start at first knot and add knots until all knots have multiplicity degree + 1
    std::vector<double>::iterator knoti = refined_knots.begin();
    while (knoti != refined_knots.end())
    {
        // Insert new knots
        int mult = degree + 1 - knot_multiplicity(*knoti);
        if (mult > 0)
        {
            std::vector<double> newKnots(mult, *knoti);
            refined_knots.insert(knoti, newKnots.begin(), newKnots.end());
        }

        // Advance to next knot
        knoti = std::upper_bound(refined_knots.begin(), refined_knots.end(), *knoti);
    }

    if (!KnotVector(refined_knots).is_regular(degree))
        throw Exception("BSplineBasis1D::refine_knots: New knot vector is not regular!");

    if (!knots.is_refinement(refined_knots))
        throw Exception("BSplineBasis1D::refine_knots: New knot vector is not a proper refinement!");

    // Return knot insertion matrix
    SparseMatrix A = build_knot_insertion_matrix(refined_knots);

    // Update knots
    knots = KnotVector(refined_knots);

    return A;
}

SparseMatrix BSplineBasis1D::build_knot_insertion_matrix(const std::vector<double> &refined_knots) const
{
    if (!KnotVector(refined_knots).is_regular(degree))
        throw Exception("BSplineBasis1D::build_knot_insertion_matrix: New knot vector is not regular!");

    if (!knots.is_refinement(refined_knots))
        throw Exception("BSplineBasis1D::build_knot_insertion_matrix: New knot vector is not a proper refinement!");

    auto n = knots.size() - degree - 1;
    auto m = refined_knots.size() - degree - 1;

    SparseMatrix A(m, n);
    //A.resize(m,n);
    A.reserve(Eigen::VectorXi::Constant(n, degree + 1));

    // Build A row-by-row
    for (unsigned int i = 0; i < m; i++)
    {
        int u = knots.index_interval(refined_knots.at(i));

        SparseMatrix R(1,1);
        R.insert(0,0) = 1;

        // For p > 0
        for (unsigned int j = 1; j <= degree; j++)
        {
            SparseMatrix Ri = build_basis_matrix(refined_knots.at(i + j), u, j);
            R = R*Ri;
        }

        // Size check
        if (R.rows() != 1 || R.cols() != (int)degree + 1)
        {
            throw Exception("BSplineBasis1D::build_knot_insertion_matrix: Incorrect matrix dimensions!");
        }

        // Insert row values
        int j = u - degree; // First insertion index
        for (int k = 0; k < R.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(R, k); it; ++it)
        {
            if (!assert_near(it.value(), .0))
                A.insert(i, j + it.col()) = it.value();
        }
    }

    A.makeCompressed();

    return A;
}

/*
 * The B-spline domain is the half-open domain [ knots.first(), knots.end() ).
 * The hack checks if x is at the right boundary (if x = knots.end()), if so,
 * a small number is subtracted from x, moving x into the half-open domain.
 */
double BSplineBasis1D::support_hack(double x) const
{
    if (x == knots.back())
        return std::nextafter(x, std::numeric_limits<double>::lowest());
    return x;
}

SparseMatrix BSplineBasis1D::reduce_support(double lb, double ub)
{
    // Check bounds
    if (lb < knots.front() || ub > knots.back())
        throw Exception("BSplineBasis1D::reduce_support: Cannot increase support!");

    unsigned int k = degree + 1;

    auto index_lower = index_supported_basis_functions(lb).front();
    auto index_upper = index_supported_basis_functions(ub).back();

    // Check lower bound index
    if (k != knot_multiplicity(knots.at(index_lower)))
    {
        int suggested_index = index_lower - 1;
        if (0 <= suggested_index)
        {
            index_lower = suggested_index;
        }
        else
        {
            throw Exception("BSplineBasis1D::reduce_support: Suggested index is negative!");
        }
    }

    // Check upper bound index
    if (knot_multiplicity(ub) == k && knots.at(index_upper) == ub)
    {
        index_upper -= k;
    }

    // New knot vector
    std::vector<double> si(knots.cbegin()+index_lower, knots.cbegin()+index_upper+k+1);

    // Construct selection matrix A
    auto num_old = knots.size()-k; // Current number of basis functions
    auto num_new = si.size()-k; // Number of basis functions after update

    if (num_old < num_new)
        throw Exception("BSplineBasis1D::reduce_support: Number of basis functions is increased instead of reduced!");

    DenseMatrix Ad = DenseMatrix::Zero(num_new, num_old);
    Ad.block(0, index_lower, num_new, num_new) = DenseMatrix::Identity(num_new, num_new);

    // Update knots
    knots = si;

    return Ad.sparseView();
}

unsigned int BSplineBasis1D::get_num_basis_functions() const
{
    return knots.size() - (degree + 1);
}

unsigned int BSplineBasis1D::get_num_basis_functions_target() const
{
    return target_num_basis_functions;
}

// Return indices of supporting basis functions at x
std::vector<unsigned int> BSplineBasis1D::index_supported_basis_functions(double x) const
{
    if (!is_supported(x))
        throw Exception("BSplineBasis1D::index_supported_basis_functions: x not inside support!");

    std::vector<unsigned int> supported;
    for (unsigned int i = 0; i < get_num_basis_functions(); ++i)
    {
        // Support of basis function i
        if (knots.at(i) <= x && x < knots.at(i+degree+1))
        {
            supported.push_back(i);

            if (supported.size() == degree + 1)
                break;
        }
    }

    // Right edge case
    if (assert_near(x, knots.back()) && knot_multiplicity(knots.back()) == degree + 1 && supported.size() < degree + 1)
    {
        auto last_basis_func = get_num_basis_functions()-1;
        if (find(supported.begin(), supported.end(), last_basis_func) == supported.end())
            supported.push_back(last_basis_func);
    }

    if (supported.empty())
        throw Exception("BSplineBasis1D::index_supported_basis_functions: No supporting basis functions");

    if (supported.size() > degree + 1)
        throw Exception("BSplineBasis1D::index_supported_basis_functions: Number of supporting basis functions larger than degree + 1!");

    return supported;
}

unsigned int BSplineBasis1D::index_longest_interval(const std::vector<double> &vec) const
{
    double longest = 0;
    double interval = 0;
    unsigned int index = 0;

    for (unsigned int i = 0; i < vec.size() - 1; i++)
    {
        interval = vec.at(i+1) - vec.at(i);
        if (longest < interval)
        {
            longest = interval;
            index = i;
        }
    }
    return index;
}

} // namespace SPLINTER

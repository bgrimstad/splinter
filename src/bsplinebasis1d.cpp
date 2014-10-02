/*
 * This file is part of the Multivariate Splines library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com)
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/


#include "include/bsplinebasis1d.h"
#include <iostream>
#include <algorithm>

namespace MultivariateSplines
{

BSplineBasis1D::BSplineBasis1D(std::vector<double> &x, unsigned int degree)
    : BSplineBasis1D(x, degree, KnotVectorType::FREE)
{
}

BSplineBasis1D::BSplineBasis1D(std::vector<double> &x, unsigned int degree, KnotVectorType knotVectorType)
    : degree(degree),
      targetNumBasisfunctions(degree+1+2) // Minimum p+1
{
    if(knotVectorType == KnotVectorType::EXPLICIT)
    {
        knots = x;
    }
    else if(knotVectorType == KnotVectorType::FREE)
    {
        knots = knotVectorFree(x);
    }
    else if(knotVectorType == KnotVectorType::REGULAR)
    {
        knots = knotVectorRegular(x);
    }
    else if(knotVectorType == KnotVectorType::EQUIDISTANT)
    {
        knots = knotVectorEquidistant(x);
    }
    else
    {
        // Assuming a regular knot vector is given
        knots = x;
    }

    assert(degree > 0);
    assert(isKnotVectorRegular());
}

SparseVector BSplineBasis1D::evaluate(double x) const
{
    SparseVector basisvalues(numBasisFunctions());

    supportHack(x);

    std::vector<int> indexSupported = indexSupportedBasisfunctions(x);

    basisvalues.reserve(indexSupported.size());

    // Iterate through the nonzero basisfunctions and store functionvalues
    for(auto it = indexSupported.begin(); it != indexSupported.end(); it++)
    {
        basisvalues.insert(*it) = deBoorCox(x, *it, degree);
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

    return basisvalues;
}

SparseVector BSplineBasis1D::evaluateDerivative(double x, int r) const
{
    // Evaluate rth derivative of basis functions at x
    // Returns vector [D^(r)B_(u-p,p)(x) ... D^(r)B_(u,p)(x)]
    // where u is the knot index and p is the degree
    int p = degree;

    // Continuity requirement
    //assert(p >= r+1);
    if(!(p >= r+1))
    {
        // Return zero-gradient
        SparseVector DB(numBasisFunctions());
        return DB;
    }

    // Check for knot multiplicity here!

    supportHack(x);

    int knotIndex = indexHalfopenInterval(x);

    // Algorithm 3.18 from Lyche and Moerken 2011
    SparseMatrix B(1,1);
    B.insert(0,0) = 1;

    for(int i = 1; i <= p-r; i++)
    {
        SparseMatrix R = buildBasisMatrix(x, knotIndex, i);
        B = B*R;
    }

    for(int i = p-r+1; i <= p; i++)
    {
        SparseMatrix DR = buildBasisMatrix(x, knotIndex, i, true);
        B = B*DR;
    }
    double factorial = std::tgamma(p+1)/std::tgamma(p-r+1);
    B = B*factorial;

    assert(B.cols() == p+1);

    // From row vector to extended column vector
    SparseVector DB(numBasisFunctions());
    DB.reserve(p+1);
    int i = knotIndex-p; // First insertion index
    for(int k = 0; k < B.outerSize(); ++k)
    for(SparseMatrix::InnerIterator it(B,k); it; ++it)
    {
        DB.insert(i+it.col()) = it.value();
    }

    return DB;
}

// Old implementation of first derivative of basis functions
DenseVector BSplineBasis1D::evaluateFirstDerivative(double x) const
{
    DenseVector values;
    values.setZero(numBasisFunctions());

    supportHack(x);

    std::vector<int> supportedBasisFunctions = indexSupportedBasisfunctions(x);

    for(int i : supportedBasisFunctions)
    {
        // Differentiate basis function
        // Equation 3.35 in Lyche & Moerken (2011)
        double b1 = deBoorCox(x, i, degree-1);
        double b2 = deBoorCox(x, i+1, degree-1);

        double t11 = knots.at(i);
        double t12 = knots.at(i+degree);
        double t21 = knots.at(i+1);
        double t22 = knots.at(i+degree+1);

        (t12 == t11) ? b1 = 0 : b1 = b1/(t12-t11);
        (t22 == t21) ? b2 = 0 : b2 = b2/(t22-t21);

        values(i) = degree*(b1 - b2);
    }

    return values;
}

// Used to evaluate basis functions - alternative to the recursive deBoorCox
SparseMatrix BSplineBasis1D::buildBasisMatrix(double x, unsigned int u, unsigned int k, bool diff) const
{
    /* Build B-spline Matrix
     * R_k in R^(k,k+1)
     * or, if diff = true, the differentiated basis matrix
     * DR_k in R^(k,k+1)
     */

    if(!(k >= 1 && k <= getBasisDegree()))
    {
        throw Exception("BSplineBasis1D::buildBasisMatrix: Incorrect input paramaters!");
    }

//    assert(u >= basisDegree + 1);
//    assert(u < ks.size() - basisDegree);

    unsigned int rows = k;
    unsigned int cols = k+1;
    SparseMatrix R(rows,cols);
    R.reserve(Eigen::VectorXi::Constant(cols,2));

    for(unsigned int i = 0; i < rows; i++)
    {
        double dk = knots.at(u+1+i) - knots.at(u+1+i-k);
        if (dk == 0)
        {
            continue;
        }
        else
        {
            if(diff)
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
                if(a != 0)
                    R.insert(i,i) = a;

                // Insert super-diagonal element
                double b = (x - knots.at(u+1+i-k))/dk;
                if(b != 0)
                    R.insert(i,i+1) = b;
            }
        }
    }

    R.makeCompressed();

    return R;
}

double BSplineBasis1D::deBoorCox(double x, int i, int k) const
{
    if(k == 0)
    {
        if(inHalfopenInterval(x, knots.at(i), knots.at(i+1)))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        double s1,s2,r1,r2;

        s1 = deBoorCoxCoeff(x, knots.at(i),   knots.at(i+k));
        s2 = deBoorCoxCoeff(x, knots.at(i+1), knots.at(i+k+1));

        r1 = deBoorCox(x, i,   k-1);
        r2 = deBoorCox(x, i+1, k-1);

        return s1*r1 + (1-s2)*r2;
    }
}

double BSplineBasis1D::deBoorCoxCoeff(double x, double x_min, double x_max) const
{
    if(x_min < x_max && x_min <= x && x <= x_max)
    {
        return (x - x_min)/(x_max - x_min);
    }
    return 0;
}

// Insert knots and compute knot insertion matrix (to update control points)
bool BSplineBasis1D::insertKnots(SparseMatrix &A, double tau, unsigned int multiplicity)
{
    if(!insideSupport(tau) || knotMultiplicity(tau) + multiplicity > degree + 1)
        return false;

    // New knot vector
    int index = indexHalfopenInterval(tau);

    std::vector<double> extKnots = knots;
    for(unsigned int i = 0; i < multiplicity; i++)
        extKnots.insert(extKnots.begin()+index+1, tau);

    assert(isKnotVectorRegular(extKnots));

    // Return knot insertion matrix
    if(!buildKnotInsertionMatrix(A, extKnots))
        return false;

    // Update knots
    knots = extKnots;

    return true;
}

bool BSplineBasis1D::refineKnots(SparseMatrix &A)
{
    // Build refine knot vector
    std::vector<double> refinedKnots = knots;

    unsigned int targetNumKnots = targetNumBasisfunctions + degree + 1;
    while(refinedKnots.size() < targetNumKnots)
    {
        int index = indexLongestInterval(refinedKnots);
        double newKnot = (refinedKnots.at(index) + refinedKnots.at(index+1))/2.0;
        refinedKnots.insert(lower_bound(refinedKnots.begin(), refinedKnots.end(), newKnot), newKnot);
    }

    assert(isKnotVectorRegular(refinedKnots) && isRefinement(refinedKnots));

    // Return knot insertion matrix
    if(!buildKnotInsertionMatrix(A, refinedKnots))
    {
        return false;
    }

    // Update knots
    knots = refinedKnots;

    return true;
}

bool BSplineBasis1D::buildKnotInsertionMatrix(SparseMatrix &A, const std::vector<double> &refinedKnots) const
{
    if (!isRefinement(refinedKnots))
    {
        throw Exception("BSplineBasis1D::buildKnotInsertionMatrix: New knot vector is not a proper refinement of the old!");
    }

    std::vector<double> knotsAug = refinedKnots;
    unsigned int n = knots.size() - degree - 1;
    unsigned int m = knotsAug.size() - degree - 1;

    //SparseMatrix A(m,n);
    A.resize(m,n);
    A.reserve(Eigen::VectorXi::Constant(n,degree+1));

    // Build A row-by-row
    for(unsigned int i = 0; i < m; i++)
    {
        int u = indexHalfopenInterval(knotsAug.at(i));

        // Assuming that p > 0
        SparseMatrix R(1,1);
        R.insert(0,0) = 1;
        for(unsigned int j = 1; j <= degree; j++)
        {
            SparseMatrix Ri = buildBasisMatrix(knotsAug.at(i+j), u, j);
            R = R*Ri;
        }

        // Size check
        if(R.rows() != 1 || R.cols() != degree+1)
        {
            throw Exception("BSplineBasis1D::buildKnotInsertionMatrix: Incorrect matrix dimensions!");
        }

        // Insert row values
        int j = u-degree; // First insertion index
        for(int k = 0; k < R.outerSize(); ++k)
        for(SparseMatrix::InnerIterator it(R,k); it; ++it)
        {
            A.insert(i,j+it.col()) = it.value();
        }
    }

    A.makeCompressed();

    return true;
}

/*
 * The B-spline domain is the half-open domain
 * [ knots.first(), knots.end() ).
 * The hack checks if x is at the right boundary
 * (x = knots.end()), and if so, subtracts a small
 * number from x, moving x into the half-open domain.
 */
void BSplineBasis1D::supportHack(double &x) const
{
    if(x == knots.back())
        x = std::nextafter(x, std::numeric_limits<double>::lowest());
}

/*
 * Finds index i such that knots.at(i) <= x < knots.at(i+1)
 * Returns false if x is outside support
 */
int BSplineBasis1D::indexHalfopenInterval(double x) const
{
    if(x < knots.front() || x > knots.back())
    {
        throw Exception("BSplineBasis1D::indexHalfopenInterval: x outside knot interval!");
    }

    // Find first knot that is larger than x
    std::vector<double>::const_iterator it = std::upper_bound(knots.begin(), knots.end(), x);

    // Return index
    int index = it - knots.begin();
    return index - 1;
}

bool BSplineBasis1D::reduceSupport(double lb, double ub, SparseMatrix &A)
{
    // Check bounds
    if(lb < knots.front() || ub > knots.back())
    {
        return false;
    }

    unsigned int k = degree + 1;

    int index_lower = indexSupportedBasisfunctions(lb).front();
    int index_upper = indexSupportedBasisfunctions(ub).back();

    // Check lower bound index
    unsigned int count = knotMultiplicity(knots.at(index_lower));
    bool is_p_regular = (k == count);

    if(!is_p_regular)
    {
        int suggested_index = index_lower - 1;
        if(0 <= suggested_index)
        {
            index_lower = suggested_index;
        }
        else
        {

#ifndef NDEBUG
            std::cout << "\n\n----------------adjust_index_for_domain_reduction-----------------" << std::endl;
            std::cout << "Error: not enough knots to guarantee controlpoint convergence" << std::endl;
            std::cout << "----------------adjust_index_for_domain_reduction-----------------\n\n" << std::endl;
#endif // NDEBUG

            return false;
        }
    }

    // Check upper bound index
    if(knotMultiplicity(ub) == k && knots.at(index_upper) == ub)
    {
        index_upper -= k;
    }

    // New knot vector
    std::vector<double> si;
    si.insert(si.begin(), knots.begin()+index_lower, knots.begin()+index_upper+k+1);

    // Construct selection matrix A
    int n_old = knots.size()-k; // Current number of basis functions
    int n_new = si.size()-k; // Number of basis functions after update

    if (n_old < n_new) return false;

    DenseMatrix Ad = DenseMatrix::Zero(n_old, n_new);
    Ad.block(index_lower, 0, n_new, n_new) = DenseMatrix::Identity(n_new, n_new);
    A = Ad.sparseView();

    // Update knots
    knots = si;

    return true;
}

double BSplineBasis1D::getKnotValue(unsigned int index) const
{
    if(index >= knots.size())
    {
        throw Exception("BSplineBasis1D:getKnotValue: Invalid knot index - Out of Range");
    }

    return knots.at(index);
}

unsigned int BSplineBasis1D::knotMultiplicity(double tau) const
{
    return std::count(knots.begin(), knots.end(), tau);
}

bool BSplineBasis1D::inHalfopenInterval(double x, double x_min, double x_max) const
{
    return (x_min <= x) && (x < x_max);
}

bool BSplineBasis1D::insideSupport(double x) const
{
    return (knots.front() <= x) && (x <= knots.back());
}

unsigned int BSplineBasis1D::numBasisFunctions() const
{
    return knots.size() - (degree + 1);
}

unsigned int BSplineBasis1D::numBasisFunctionsTarget() const
{
    return targetNumBasisfunctions;
}

// Return indices of supporting basis functions at x
std::vector<int> BSplineBasis1D::indexSupportedBasisfunctions(double x) const
{
    std::vector<int> ret;
    if(insideSupport(x))
    {
        int last = indexHalfopenInterval(x);
        if(last < 0)
        {
            last = knots.size() - 1 - (degree + 1);
        }
        int first = std::max((int)(last - degree), 0);
        for(int i = first; i <= last; i++)
        {
            ret.push_back(i);
        }
    }
    return ret;
}

unsigned int BSplineBasis1D::indexLongestInterval() const
{
    return indexLongestInterval(knots);
}

unsigned int BSplineBasis1D::indexLongestInterval(const std::vector<double> &vec) const
{
    double longest = 0;
    double interval = 0;
    unsigned int index = 0;

    for(unsigned int i = 0; i < vec.size() - 1; i++)
    {
        interval = vec.at(i+1) - vec.at(i);
        if(longest < interval)
        {
            longest = interval;
            index = i;
        }
    }
    return index;
}

bool BSplineBasis1D::isKnotVectorRegular() const
{
    return isKnotVectorRegular(knots);
}

bool BSplineBasis1D::isKnotVectorRegular(const std::vector<double> &vec) const
{
    // Check size
    if(vec.size() < 2*(degree+1))
        return false;

    // Check first knots
    if(std::count(vec.begin(), vec.begin()+degree+1, vec.front()) != degree+1)
        return false;

    // Check last knots
    if(std::count(vec.end()-degree-1, vec.end(), vec.back()) != degree+1)
        return false;

    // Check order
    if(!std::is_sorted(vec.begin(), vec.end()))
        return false;

    // Check multiplicity of knots
    for(std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
    {
        if(count(vec.begin(), vec.end(), *it) > degree+1)
            return false;
    }

    return true;
}

bool BSplineBasis1D::isRefinement(const std::vector<double> &refinedKnots) const
{
    // Check size
    if(refinedKnots.size() < knots.size())
        return false;

    // Check that t is regular
    if(!isKnotVectorRegular(refinedKnots))
        return false;

    // Check that each element in knots occurs at least as many times in refinedKnots
    for(std::vector<double>::const_iterator it = knots.begin() ; it != knots.end(); ++it)
    {
        int m_tau = count(knots.begin(), knots.end(), *it);
        int m_t = count(refinedKnots.begin(), refinedKnots.end(), *it);
        if (m_t < m_tau) return false;
    }

    // Check that range is not changed
    if(knots.front() != refinedKnots.front()) return false;
    if(knots.back() != refinedKnots.back()) return false;

    return true;
}

/*
 * Creates a regular knot vector of n+p+1 equidistant knots,
 * where the first and last knot is repeated p+1 times,
 * and n is the number of unique values in x.
 *
 * This knot vector can be used for all degrees > 0.
 */
std::vector<double> BSplineBasis1D::knotVectorEquidistant(std::vector<double> &X) const
{
    // Copy X -> sort -> remove duplicates -> resize = a sorted vector of unique values
    std::vector<double> uniqueX(X);
    sort(uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy(uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize(distance(uniqueX.begin(),it));

    // Minimum number of interpolation points
    if(uniqueX.size() < degree+1)
    {
        std::ostringstream e;
        e << "BSplineBasis1D::knotVectorFree: Only " << uniqueX.size()
          << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
          << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    std::vector<double> knots;
    for(unsigned int i = 0; i < degree; i++)
        knots.push_back(uniqueX.front());

    std::vector<double> equiKnots = linspace(uniqueX.front(), uniqueX.back(), uniqueX.size() - degree + 1);
    for(auto it = equiKnots.begin(); it != equiKnots.end(); ++it)
        knots.push_back(*it);

    for(unsigned int i = 0; i < degree; i++)
        knots.push_back(uniqueX.back());

    return knots;
}

/*
 * Repeats first and last knot p+1 times. Removes no knots.
 */
std::vector<double> BSplineBasis1D::knotVectorRegular(std::vector<double> &X) const
{
    // Copy X -> sort -> remove duplicates -> resize = a sorted vector of unique values
    std::vector<double> uniqueX(X);
    sort(uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy(uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize(distance(uniqueX.begin(),it));

    std::vector<double> knots;
    it = uniqueX.begin();

    // Repeat first knot p + 1 times (for interpolation of start point)
    for(unsigned int i = 0; i < degree + 1; i++)
    {
        knots.push_back(*it);
    }

    // Add unique knots
    for(it = uniqueX.begin()+1; it < uniqueX.end()-1; it++)
    {
        knots.push_back(*it);
    }

    // Repeat last knot p + 1 times (for interpolation of end point)
    it = uniqueX.end()-1; // Last element in uniqueX
    for(unsigned int i = 0; i < degree + 1; i++)
    {
        knots.push_back(*it);
    }

    // Number of knots in a (p+1)-regular knot vector
    assert(knots.size() == uniqueX.size() + 2*degree);

    return knots;
}

/*
 * Free knot vector for degrees 1, 2, and 3 only.
 * The free knot vector ensure that the first and last knot is repeated p + 1 times,
 * and that the total number of knots is n + p + 1.
 * Example for degree=3 and x={a,b,c,d,e,f,g,h}: knots={a,a,a,a,c,d,e,f,h,h,h,h}
 */
std::vector<double> BSplineBasis1D::knotVectorFree(std::vector<double> &X) const
{
    // Copy X -> sort -> remove duplicates -> resize = a sorted vector of unique values
    std::vector<double> uniqueX(X);
    sort(uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy(uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize(distance(uniqueX.begin(),it));

    // The minimum number of samples from which a free knot vector can be created
    if(uniqueX.size() < degree+1)
    {
        std::ostringstream e;
        e << "BSplineBasis1D::knotVectorFree: Only " << uniqueX.size()
          << " unique interpolation points are given. A minimum of degree+1 = " << degree+1
          << " unique points are required to build a B-spline basis of degree " << degree << ".";
        throw Exception(e.str());
    }

    std::vector<double> knots;
    it = uniqueX.begin();

    // Repeat first x value p + 1 times (for interpolation of start point)
    for(unsigned int i = 0; i < degree + 1; i++)
    {
        knots.push_back(*it);
    }

    if(degree == 1)
    {
        // No knots removed
        for (it = uniqueX.begin()+1; it < uniqueX.end()-1; it++)
        {
            knots.push_back(*it);
        }
    }
    else if(degree == 2)
    {
        // First knot removed
        for (it = uniqueX.begin()+2; it < uniqueX.end()-1; it++)
        {
            knots.push_back(*it);
        }
    }
    else if(degree == 3)
    {
        // First and last knot removed
        for (it = uniqueX.begin()+2; it < uniqueX.end()-2; it++)
        {
            knots.push_back(*it);
        }
    }
    else
    {
        throw Exception("BSplineBasis1D::knotVectorFree: A free knot vector is only supported by degrees 1, 2, and 3!");
    }

    // Repeat last x value p + 1 times (for interpolation of end point)
    it = uniqueX.end()-1; // Last element in uniqueX
    for(unsigned int i = 0; i < degree + 1; i++)
    {
        knots.push_back(*it);
    }

    // Number of knots in a (p+1)-regular knot vector
    // Ensures a square matrix when calculating the control points
    assert(knots.size() == uniqueX.size() + degree + 1);

    return knots;
}

std::vector<double> BSplineBasis1D::linspace(double start, double stop, unsigned int points) const
{
    std::vector<double> ret;
    double dx = 0;
    if(points > 1)
        dx = (stop - start)/(points-1);
    for(unsigned int i = 0; i < points; ++i)
        ret.push_back(start + i*dx);
    return ret;
}

} // namespace MultivariateSplines

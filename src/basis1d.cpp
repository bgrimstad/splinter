#include "basis1d.h"
#include <iostream>
#include <algorithm>
//#include "timer.h"

//using namespace std;
using std::cout;
using std::endl;

Basis1D::Basis1D(std::vector<double> &x, int degree)
    : Basis1D(x, degree, KnotSequenceType::FREE)
{
}

Basis1D::Basis1D(std::vector<double> &x, int degree, KnotSequenceType knotSequenceType)
    : degree(degree),
      targetNumBasisfunctions(degree+1+2) // Minimum p+1
{
    if (knotSequenceType == KnotSequenceType::EXPLICIT)
    {
        knots = x;
    }
    else if (knotSequenceType == KnotSequenceType::FREE)
    {
        knots = knotSequenceFree(x);
    }
    else if (knotSequenceType == KnotSequenceType::REGULAR)
    {
        knots = knotSequenceRegular(x);
    }
    else
    {
        // Assuming a regular knot sequence is given
        knots = x;
    }

    assert(degree > 0);
    assert(isKnotSequenceRegular());
}

SparseVector Basis1D::evaluate(const double x) const
{
    SparseVector basisvalues(numBasisFunctions());

    //the following assumes regular sequence: test if x is at the end of the support,
    //because the only the first/last basis is nonzero (and equal to one)
    if (knots.front() == x)
    {
        // if evaluated at the end of the last interval
        // set the last basisfunction equal to 1, all others are 0
        basisvalues.reserve(1);
        basisvalues.insert(0) = 1;
        return basisvalues;
    }

    if (knots.back() == x)
    {
        // if evaluated at the end of the last interval
        // set the last basisfunction equal to 1, all others are 0
        basisvalues.reserve(1);
        basisvalues.insert(numBasisFunctions() - 1) = 1;
        return basisvalues;
    }

    std::vector<int> indexSupported = indexSupportedBasisfunctions(x);

    basisvalues.reserve( indexSupported.size() );

    // Iterate through the nonzero basisfunctions and store functionvalues
    for (auto itr = indexSupported.begin(); itr != indexSupported.end(); itr++)
    {
        basisvalues.insert(*itr) = deBoorCox(x, *itr, degree);
    }

//    // Must have regular knot sequence to do this - alternative evaluation using basis matrix
//    // Do knot hack
//    int knotIndex = indexHalfopenInterval(x); // knot index
//    if (knotIndex > 0)
//    {
//        SparseMatrix basisvalues2 = buildBsplineMatrix(x, knotIndex, 1);
//        for (int i = 2; i <= basisDegree; i++)
//        {
//            SparseMatrix Ri = buildBsplineMatrix(x, knotIndex, i);
//            basisvalues2 = basisvalues2*Ri;
//        }
//        basisvalues2.makeCompressed();

//        assert(basisvalues2.rows() == 1);
//        assert(basisvalues2.cols() == basisDegree + 1);

//    }
//    else
//    {
//        cout << "Knot index not good!" << endl;
//        cout << knotIndex << endl;
//        cout << x << " knotseq " << knots.front() << "/" << knots.back() << endl;
//        exit(1);
//    }

    return basisvalues;
}

SparseVector Basis1D::evaluateDerivative(double x, int r) const
{
    // Evaluate rth derivative of basis functions at x
    // Returns vector [D^(r)B_(u-p,p)(x) ... D^(r)B_(u,p)(x)]
    // where u is the knot index and p is the degree
    int p = degree;

    // Continuity requirement
    //assert(p >= r+1);
    if (!(p >= r+1))
    {
        // Return zero-gradient
        SparseVector DB(numBasisFunctions());
        return DB;
    }

    // Check for knot multiplicity here!

    // Knot hack
    if (knots.back() == x) x = x-1e-12;

    int knotIndex = indexHalfopenInterval(x);

    if (0 <= knotIndex) // -1 if the x value was outside the knot range
    {
        // Algorithm 3.18 from Lyche and Moerken 2011
        SparseMatrix B(1,1);
        B.insert(0,0) = 1;

        for (int i = 1; i <= p-r; i++)
        {
            SparseMatrix R = buildBasisMatrix(x, knotIndex, i);
            B = B*R;
        }

        for (int i = p-r+1; i <= p; i++)
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
        int i = knotIndex-p;
        for (int k = 0; k < B.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(B,k); it; ++it)
        {
            if (it.value() != 0)
                DB.insert(i) = it.value();
            i++;
        }

        return DB;
    }
    else
    {
        cout << "Knot index not good!" << endl;
        cout << knotIndex << endl;
        cout << x << " knotseq " << knots.front() << "/" << knots.back() << endl;
        cout << "Delta knot: " << std::setprecision(10) << knots.back()-knots.front() << endl;
        SparseVector DB(numBasisFunctions());
        return DB;
    }

}

// Old implementation of first derivative of basis functions
DenseVector Basis1D::evaluateFirstDerivative(double x) const
{
    DenseVector values; values.setZero(numBasisFunctions());

    // Knot hack
    if (knots.back() == x) x = x-1e-12;

    int first_knot =  indexHalfopenInterval(x);

    if ( 0 <= first_knot ) //  -1 if the x value was outside the knot range
    {
        std::vector<int> supportedBasisFunctions = indexSupportedBasisfunctions(x);

        for (int i = supportedBasisFunctions.front(); i <= supportedBasisFunctions.back(); i++)
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
    }

    return values;
}

// Used to evaluate basis functions - alternative to the recursive deBoorCox
SparseMatrix Basis1D::buildBasisMatrix(double x, int u, int k, bool diff) const
{
    /* Build B-spline Matrix
     * R_k in R^(k,k+1)
     * or, if diff = true, the differentiated basis matrix
     * DR_k in R^(k,k+1)
     */
    assert(1 <= k);
    assert(k <= getBasisDegree());
//    assert(u >= basisDegree + 1);
//    assert(u < ks.size() - basisDegree);

    int rows = k;
    int cols = k+1;
    SparseMatrix R(rows,cols);
    R.reserve(Eigen::VectorXi::Constant(cols,2));

    for (int i = 0; i < rows; i++)
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
                //if (a != 0)
                    R.insert(i,i) = a;

                // Insert super-diagonal element
                double b = (x - knots.at(u+1+i-k))/dk;
                //if (b != 0)
                    R.insert(i,i+1) = b;
            }
        }
    }

    R.makeCompressed();

    return R;
}

double Basis1D::deBoorCox(double x, int i, int k) const
{
    if ( 0 == k )
    {
        if (inHalfopenInterval(x, knots.at(i), knots.at(i+1)))
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

double Basis1D::deBoorCoxCoeff(double x, double x_min, double x_max) const
{
    if (x_min < x_max && x_min <= x && x <= x_max)
    {
        return (x - x_min)/(x_max - x_min);
    }
    return 0;
}

// Insert knots and compute knot insertion matrix (to update control points)
bool Basis1D::insertKnots(SparseMatrix &A, double tau, unsigned int multiplicity)
{
    if (!insideSupport(tau) || knotMultiplicity(tau) + multiplicity > degree + 1)
        return false;

    // New knot vector
    int index = indexHalfopenInterval(tau);

    std::vector<double> extKnots = knots;
    for (unsigned int i = 0; i < multiplicity; i++)
        extKnots.insert(extKnots.begin()+index+1, tau);

    assert(isKnotSequenceRegular(extKnots));

    // Return knot insertion matrix
    if (!buildKnotInsertionMatrix(A, extKnots))
        return false;

    // Update knots
    knots = extKnots;

    return true;
}

bool Basis1D::refineKnots(SparseMatrix &A)
{
    // Build refine knot vector
    std::vector<double> refinedKnots = knots;

    unsigned int targetNumKnots = targetNumBasisfunctions + degree + 1;
    while (refinedKnots.size() < targetNumKnots)
    {
        int index = indexLongestInterval(refinedKnots);
        double newKnot = (refinedKnots.at(index) + refinedKnots.at(index+1))/2.0;
        refinedKnots.insert(lower_bound(refinedKnots.begin(), refinedKnots.end(), newKnot), newKnot);
    }

    assert(isKnotSequenceRegular(refinedKnots));

    // Return knot insertion matrix
    if (!buildKnotInsertionMatrix(A, refinedKnots))
        return false;

    // Update knots
    knots = refinedKnots;

    return true;
}

bool Basis1D::buildKnotInsertionMatrix(SparseMatrix &A, const std::vector<double> &refinedKnots) const
{
    if (!isRefinement(refinedKnots))
        return false;

    std::vector<double> knotsAug = refinedKnots;
    int n = knots.size() - degree - 1;
    int m = knotsAug.size() - degree - 1;

    //SparseMatrix A(m,n);
    A.resize(m,n);
    A.reserve(Eigen::VectorXi::Constant(n,degree+1));

    // Build A row-by-row
    for (int i = 0; i < m; i++)
    {
        int u = indexHalfopenInterval(knotsAug.at(i));

        // Assuming that p > 0
        SparseMatrix R(1,1);
        R.insert(0,0) = 1;
        for (int j = 1; j <= degree; j++)
        {
            SparseMatrix Ri = buildBasisMatrix(knotsAug.at(i+j), u, j);
            R = R*Ri;
        }

        // Size check
        assert(R.rows() == 1);
        assert(R.cols() == degree+1);

        // Insert row values
        int j = u-degree;
        for (int k = 0; k < R.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(R,k); it; ++it)
        {
            if (it.value() != 0)
                A.insert(i,j) = it.value(); //insert
            j++;
        }
    }

    A.makeCompressed();

    return true;
}

int Basis1D::indexHalfopenInterval(double x) const
{
    // Returns index i such that knots.at(i) <= x < knots.at(i+1)
    //returns -99 if the x value is outside the knot range.
    //returns -1  if the x value is identical to the last knot
    //returns index of first nonzero interval otherwise

    if (0 == knots.size() || x < knots.front() || x > knots.back())
    {
        return -99;
    }
    else
    {
        // Find first knot that is larger than x
        std::vector<double>::const_iterator iter = std::upper_bound(knots.begin(), knots.end(), x);
        int index = iter - knots.begin();

        if (knots.end() ==  iter)
        {
            // NOTE: x is equal to or larger than the last knot.
            //  cout << "Value equal to the last knot. Treat as special case? Returning -1" << endl;
            return -1;
        }
        else
        {
            return index - 1;
        }
    }

    return -55; // Unreachable; included to avoid warning
}

bool Basis1D::reduceSupport(double lb, double ub, SparseMatrix &A)
{
    // Check bounds
    if (lb < knots.front() || ub > knots.back())
    {
        cout << "OUCH! Outside knots!" << endl;
        return false;
    }

    int k = degree + 1;

    int index_lower = indexSupportedBasisfunctions(lb).front();
    int index_upper = indexSupportedBasisfunctions(ub).back();

    // Check lower bound index
    int count = knotMultiplicity( knots.at(index_lower) );
    bool is_p_regular = (k == count);

    if (!is_p_regular)
    {
        int suggested_index = index_lower - 1;
        if (0 <= suggested_index)
        {
            index_lower = suggested_index;
        }
        else
        {
            cout << "\n\n----------------adjust_index_for_domain_reduction-----------------" << endl;
            cout << "Error: not enough knots to guarantee controlpoint convergence" << endl;
            cout << "----------------adjust_index_for_domain_reduction-----------------\n\n" << endl;
            return false;
        }
    }

    // Check upper bound index
    if (knotMultiplicity(ub) == k && knots.at(index_upper) == ub)
    {
        index_upper -= k;
    }

    // New knot vector
    std::vector<double> si;
    si.insert(si.begin(), knots.begin()+index_lower, knots.begin()+index_upper+k+1);

    // Construct selection matrix
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

double Basis1D::getKnotValue(int index) const
{
    if (index >= (int)knots.size())
        return knots.back();

    if (index < 0)
        return knots.front();

    return knots.at(index);
}

int Basis1D::knotMultiplicity(const double& tau) const
{
    return std::count(knots.begin(), knots.end(), tau);
}

bool Basis1D::inHalfopenInterval(double x, double x_min, double x_max) const
{
    return (x_min <= x) && (x < x_max);
}

bool Basis1D::insideSupport(const double &x) const
{
    return (knots.front() <= x) && (x <= knots.back());
}

int Basis1D::numBasisFunctions() const
{
    return knots.size() - (degree + 1);
}

int Basis1D::numBasisFunctionsTarget() const
{
    return targetNumBasisfunctions;
}

// Return indices of supporting basis functions at x
std::vector<int> Basis1D::indexSupportedBasisfunctions(double x) const
{
    std::vector<int> ret;
    if(insideSupport(x))
    {
        int last = indexHalfopenInterval(x);
        if (last < 0)
        {
            last = knots.size() - 1 - (degree + 1);
        }
        int first = std::max(last - degree, 0);
        for (int i = first; i <= last; i++)
        {
            ret.push_back(i);
        }
    }
    return ret;
}

int Basis1D::indexLongestInterval() const
{
    return indexLongestInterval(knots);
}

int Basis1D::indexLongestInterval(const std::vector<double> &vec) const
{
    double longest = 0;
    double interval = 0;
    int index = 0;

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

bool Basis1D::isKnotSequenceRegular() const
{
    return isKnotSequenceRegular(knots);
}

bool Basis1D::isKnotSequenceRegular(const std::vector<double> &vec) const
{
    // Check size
    if (vec.size() < 2*(degree+1))
        return false;

    // Check first knots
    if (std::count(vec.begin(), vec.begin()+degree+1, vec.front()) != degree+1)
        return false;

    // Check last knots
    if (std::count(vec.end()-degree-1, vec.end(), vec.back()) != degree+1)
        return false;

    // Check order
    if (!std::is_sorted(vec.begin(), vec.end()))
        return false;

    // Check multiplicity of knots
    for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it)
    {
        if (count(vec.begin(), vec.end(), *it) > degree+1)
            return false;
    }

    return true;
}

bool Basis1D::isRefinement(const std::vector<double> &refinedKnots) const
{
    // Check size
    if (refinedKnots.size() < knots.size())
        return false;

    // Check that t is regular
    if (!isKnotSequenceRegular(refinedKnots))
        return false;

    // Check that each elements in tau occurs at least as man times in t
    for (std::vector<double>::const_iterator it = knots.begin() ; it != knots.end(); ++it)
    {
        int m_tau = count(knots.begin(), knots.end(), *it);
        int m_t = count(refinedKnots.begin(), refinedKnots.end(), *it);
        if (m_t < m_tau) return false;
    }

    // Check that range is not changed
    if (knots.front() != refinedKnots.front()) return false;
    if (knots.back() != refinedKnots.back()) return false;

    return true;
}

std::vector<double> Basis1D::knotSequenceRegular(std::vector<double> &X)
{
    // Copy X -> sort -> remove duplicates -> resize = a sorted vector of unique values
    std::vector<double> uniqueX(X);
    sort (uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy (uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize( distance(uniqueX.begin(),it) );

    // Multiplicity of the first and last knot.
    // basisDegree + 1 creates discontinuity that allows for interpolation at endpoints
    int repeats_at_ends = degree + 1;

    // Number of knots in a (p+1)-regular knot sequence
    int num_knots = uniqueX.size() + 2*degree;

    std::vector<double> knots;
    it = uniqueX.begin();

    // Repeat first knot p + 1 times
    for (int i = 0; i < repeats_at_ends; i++)
    {
        knots.push_back(*it);
    }

    // Add unique knots
    for (it = uniqueX.begin()+1; it < uniqueX.end()-1; it++)
    {
        knots.push_back(*it);
    }

    // Repeat last knot p + 1 times
    it = uniqueX.end()-1; // Last element in uniqueX
    for (int i = 0; i < repeats_at_ends; i++)
    {
        knots.push_back(*it);
    }

    assert(knots.size() == (unsigned)num_knots);

    return knots;
}

// Only implemented for cubic splines!
std::vector<double> Basis1D::knotSequenceFree(std::vector<double> &X)
{
    // Copy X -> sort -> remove duplicates -> resize = a sorted vector of unique values
    std::vector<double> uniqueX(X);
    sort (uniqueX.begin(), uniqueX.end());
    std::vector<double>::iterator it = unique_copy (uniqueX.begin(), uniqueX.end(), uniqueX.begin());
    uniqueX.resize( distance(uniqueX.begin(),it) );

    // Multiplicity of the first and last knot.
    // basisDegree + 1 creates discontinuity that allows for interpolation at endpoints
    int repeats_at_ends = degree + 1;

    // Number of knots in a (p+1)-regular knot sequence
    int num_knots = uniqueX.size() + repeats_at_ends;

    // NOTE: when splineDegree = 3 and uniqueX.size() = 2 (two samples)
    // num_knots = 6, although 8 knots are required for a p-regular knot sequence!
    // The code fail in this special case.
    // if (splineDegree == 3 && num_knots < 8) num_knots = 8; // Does not work

    std::vector<double> knots;
    it = uniqueX.begin();

    // Repeat first x value p + 1 times, skip second first and second last, repeat last p + 1 times
    for (int i = 0; i < repeats_at_ends; i++)
    {
        knots.push_back(*it);
    }

    if (1 == degree)
    {
        knots.insert(knots.end(), uniqueX.begin()+1, uniqueX.end()-1);
    }
    else if( 3 == degree)
    {
        for (it = uniqueX.begin()+2 ; it < uniqueX.end()-2; it++)
        {
            knots.push_back(*it);
        }
    }
    else
    {
        cout << "Only 1. and 3. degree supported for knot generation." << endl;
    }

    // Repeat knots at end
    it = uniqueX.end()-1; // Last element in uniqueX
    for (int i = 0; i < repeats_at_ends; i++)
    {
        knots.push_back(*it);
    }

    assert(knots.size() == (unsigned)num_knots);

    return knots;
}

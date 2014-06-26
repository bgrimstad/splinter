#include "bspline.h"
#include "tensorindex.h"
#include "unsupported/Eigen/KroneckerProduct"
#include "Eigen/IterativeLinearSolvers"
#include "basis.h"
#include "basis1d.h"
//#include "timer.h"

#include <iostream>

using std::cout;
using std::endl;

// Constructor for multivariate B-splines (with one or more outputs)
Bspline::Bspline(DenseMatrix coefficients, std::vector< std::vector<double> > knotSequences, std::vector<int> basisDegrees)
    : coefficients(coefficients)
{
    numInputs = knotSequences.size();
    numOutputs = coefficients.rows();

    basis = Basis(knotSequences, basisDegrees, KnotSequenceType::EXPLICIT);

    calculateKnotAverages();

    init();

    checkControlPoints();
}

// Constructors for interpolation data (InterpolationTable)
//Bspline::Bspline(InterpolationTable &data, int basisdegree)
//    : Bspline(data, std::vector<int>(data.getDimX(), basisdegree ))
//{
//    // All bases have the same degree
//}

//Bspline::Bspline(InterpolationTable &data, std::vector<int> basisdegree)
//{
//    // Check data
//    assert(data.isGridComplete());

//    std::vector< std::vector<double> > xdata = data.getTransposedTableX();
//    std::vector< std::vector<double> > ydata = data.getTransposedTableY();

//    numInputs = xdata.size();
//    numOutputs = ydata.size();

//    // Construct basis and control points
//    //basis = Basis(xdata, basisdegree);
//    //calculateControlPoints(xdata, ydata);

//    BsplineType bsplineType = BsplineType::FREE;

//    if (bsplineType == BsplineType::PSPLINE)
//    {
//        basis = Basis(xdata, basisdegree, KnotSequenceType::FREE);
//        calculateControlPointsPspline(xdata, ydata);
//    }
//    else
//    {
//        // Default is FREE
//        basis = Basis(xdata, basisdegree, KnotSequenceType::FREE);
//        calculateControlPoints(xdata, ydata);
//    }

//    init();

//    checkControlPoints();
//}

// Constructors for interpolation data (SortedDataTable)
Bspline::Bspline(SortedDataTable &data, int basisdegree)
    : Bspline(data, std::vector<int>(data.getXDimension(), basisdegree ))
{
    // All bases have the same degree
}

Bspline::Bspline(SortedDataTable &data, std::vector<int> basisdegree)
{
    // Check data
    assert(data.isGridComplete());

    std::vector< std::vector<double> > xdata = data.getTransposedTableX();
    std::vector< std::vector<double> > ydata = data.getTransposedTableY();

    numInputs = xdata.size();
    numOutputs = ydata.size();

    // Construct basis and control points
    //basis = Basis(xdata, basisdegree);
    //calculateControlPoints(xdata, ydata);

    BsplineType bsplineType = BsplineType::FREE;

    if (bsplineType == BsplineType::PSPLINE)
    {
        basis = Basis(xdata, basisdegree, KnotSequenceType::FREE);
        calculateControlPointsPspline(xdata, ydata);
    }
    else
    {
        // Default is FREE
        basis = Basis(xdata, basisdegree, KnotSequenceType::FREE);
        calculateControlPoints(xdata, ydata);
    }

    init();

    checkControlPoints();
}

void Bspline::init()
{
    bool initialKnotRefinement = false;
    if (initialKnotRefinement)
    {
        refineKnotSequences();
        checkControlPoints();
    }
}

DenseVector Bspline::evaluate(DenseVector &x)
{
    if (!valueInsideDomain(x))
    {
        DenseVector y(numOutputs);
        y.setZero(numOutputs);
        //TODO return
    }

    SparseVector tensorvalues = basis.evaluate(x);
    return coefficients*tensorvalues;
}

DenseMatrix Bspline::jacobian(DenseVector &x)
{
    if (valueInsideDomain(x))
    {
//        Timer timer;
//        timer.start();

        DenseMatrix Bi = basis.evaluateBasisJacobianOld(x);
//        timer.stop();
//        cout << "Time old Jacobian: " << timer.getMicroSeconds() << endl;
//        timer.reset();
//        timer.start();
//        SparseMatrix Ji = basis.evaluateBasisJacobian(x);
//        timer.stop();
//        cout << "Time new Jacobian: " << timer.getMicroSeconds() << endl;
        // New Jacobian is about 5 times slower at this point...

//        // Test difference in Jacobians
//        DenseMatrix dJ = Bi - Ji;
//        DenseVector errorVec = dJ.rowwise().maxCoeff();
//        DenseVector error = errorVec.colwise().maxCoeff();

//        if (abs(error(0)) > 1e-10) cout << "NOTABLE DIFFERENCE IN JACOBIANS: " << abs(error(0)) << endl;
//        cout << abs(error(0)) << endl;
//        assert(abs(error(0)) < 1e-10);

        return coefficients*Bi;
    }
    else
    {
        cout << "Warning: returning empty Jacobian!" << endl;
        exit(1);
        DenseMatrix y; y.setZero(numOutputs, numInputs);
        return y;
    }
}

DenseMatrix Bspline::hessian(DenseVector &x)
{
    DenseMatrix H; H.setZero(1,1);
    if (numOutputs == 1)
    {
        // TODO: use SparseMatrix and implement for numOutputs > 1
        DenseMatrix identity = DenseMatrix::Identity(numInputs,numInputs);
//        cout << "Identity matrix: " << endl;
//        cout << identity.rows() << "/" << identity.cols() << endl;
        DenseMatrix caug;
//        cout << "Control point vector: " << endl;
//        cout << coefficients.rows() << "/" << coefficients.cols() << endl;
        caug = kroneckerProduct(identity, coefficients);
//        cout << "Augmented control point vector: " << endl;
//        cout << caug.rows() << "/" << caug.cols() << endl;
        DenseMatrix DB = basis.evaluateBasisHessian(x);
//        cout << "Basis Hessian: " << endl;
//        cout << DB.rows() << "/" << DB.cols() << endl;
        H = caug*DB;
        return H;
    }
    else
    {
        cout << "Warning: returning empty Hessian!" << endl;
        return H;
    }
}

std::vector< std::vector<double> > Bspline::getKnotVectors() const
{
    return basis.getKnotVectors();
}

std::vector<double> Bspline::getDomainUpperBound() const
{
    return basis.getSupportUpperBound();
}

std::vector<double> Bspline::getDomainLowerBound() const
{
    return basis.getSupportLowerBound();
}

DenseMatrix Bspline::getControlPoints() const
{
    int nc = coefficients.cols();
    DenseMatrix controlPoints(numInputs + numOutputs, nc);

    controlPoints.block(0,          0, numInputs,   nc) = knotaverages;
    controlPoints.block(numInputs,  0, numOutputs,  nc) = coefficients;

    return controlPoints;
}

void Bspline::setControlPoints(DenseMatrix &controlPoints)
{
    assert(controlPoints.rows() == numInputs + numOutputs);
    int nc = controlPoints.cols();

    knotaverages = controlPoints.block(0,          0, numInputs,   nc);
    coefficients = controlPoints.block(numInputs,  0, numOutputs,  nc);

    checkControlPoints();
}

bool Bspline::checkControlPoints() const
{
    assert(coefficients.cols() == knotaverages.cols());
    assert(knotaverages.rows() == numInputs);
    assert(coefficients.rows() == numOutputs);
    return true;
}

bool Bspline::valueInsideDomain(DenseVector x)
{
    return basis.valueInsideSupport(x);
}

bool Bspline::reduceDomain(std::vector<double> lb, std::vector<double> ub, bool regularKnotsequences, bool refineKnotsequences)
{
    if (lb.size() != numInputs || ub.size() != numInputs)
        return false;

    std::vector<double> sl = basis.getSupportLowerBound();
    std::vector<double> su = basis.getSupportUpperBound();

    bool isStrictSubset = false;
    double minDistance = 0; //1e-6; // TODO: Set to zero, B-spline constraint class is responsible

    for (unsigned int dim = 0; dim < numInputs; dim++)
    {
        if (ub.at(dim) - lb.at(dim) <  minDistance)
        {
            cout << "Cannot reduce B-spline domain: inconsistent bounds!" << endl;
            return false;
        }
        if (su.at(dim) - ub.at(dim) > minDistance)
        {
            isStrictSubset = true;
            su.at(dim) = ub.at(dim);
        }

        if (lb.at(dim) - sl.at(dim) > minDistance)
        {
            isStrictSubset = true;
            sl.at(dim) = lb.at(dim);
        }
    }

    if (isStrictSubset)
    {
        if (regularKnotsequences)
        {
            if (!regularSequences(sl, su))
            {
                cout << "Failed to regularize knot vectors!" << endl;
                exit(1);
                return false;
            }
        }

        // Remove knots and control points that are unsupported with the new bounds
        if (!removeUnsupportedBasisFunctions(sl, su))
        {
            cout << "Failed to remove unsupported basis functions!" << endl;
            exit(1);
        }

        // Refine knots
        if (refineKnotsequences)
        {
            if (!refineKnotSequences())
            {
                cout << "Failed to refine knot vectors!" << endl;
                exit(1);
                return false;
            }
        }
    }

    return true;
}

// Calculates knot averages: assumes that basis is initialized!
void Bspline::calculateKnotAverages()
{
    // Calculate knot averages for each knot sequence
    std::vector< std::vector<double> > knotAverages;
    for (unsigned int i = 0; i < numInputs; i++)
    {
        std::vector<double> knots = basis.getKnotVector(i);
        std::vector<double> knotAvgs(basis.numBasisFunctions(i));
        for (int j = 0; j < basis.numBasisFunctions(i); j++)
        {
            double knotAvg = 0;
            for (int k = j+1; k <= j+basis.getBasisDegree(i); k++)
            {
                knotAvg += knots.at(k);
            }
            knotAvgs.at(j) = knotAvg/basis.getBasisDegree(i);
        }
        knotAverages.push_back(knotAvgs);
    }

    // Fill out knot average matrix
    TensorIndex ti(basis.getTensorIndexDimension());
    assert(ti.vectorSize() == basis.numBasisFunctions());
    knotaverages.resize(numInputs, ti.vectorSize());

    for (int i = 0; i < ti.vectorSize(); i++)
    {
        std::vector<int> idx = ti.vectorToTensorIndex(i);

        for (unsigned int j = 0; j < idx.size(); j++)
        {
            knotaverages(j,i) = knotAverages.at(j).at(idx.at(j));
        }
    }
}

void Bspline::calculateControlPoints(std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y)
{
    // Input data tests
    assert(X.size() > 0);
    assert(Y.size() > 0);

    // Assuming regular grid
    unsigned int numSamples = X.at(0).size();

    for (unsigned int i = 0; i < numInputs; i++)
    {
        assert(X.at(i).size() == numSamples);
    }

    for (unsigned int i = 0; i < numOutputs; i++)
    {
        assert(Y.at(i).size() == numSamples);
    }

    /* Setup and solve equations Ac = b,
     * A = basis functions at sample x-values,
     * b = sample y-values when calculating control coefficients,
     * b = sample x-values when calculating knot averages
     * c = control coefficients or knot averages.
     */

    SparseMatrix A;
    DenseMatrix Bx, By;
    controlPointEquationLHS(X, A);
    controlPointEquationRHS(Y, By);
    controlPointEquationRHS(X, Bx);

    DenseMatrix CY;
    DenseMatrix CX;

    int numEquations = A.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if (!solveAsDense)
    {
        cout << "Computing B-spline control points using sparse solver." << endl;
        bool successfulSolve = (solveSparseLSE(A,Bx,CX) && solveSparseLSE(A,By,CY));

        solveAsDense = !successfulSolve;
    }

    if (solveAsDense)
    {
        cout << "Computing B-spline control points using dense solver." << endl;
        DenseMatrix Ad = A.toDense();
        bool successfulSolve = (solveDenseLSE(Ad,Bx,CX) && solveDenseLSE(Ad,By,CY));
        assert(successfulSolve);
    }

    coefficients = CY.transpose();
    knotaverages = CX.transpose();
}

void Bspline::calculateControlPointsPspline(std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y)
{
    // Input data tests
    assert(X.size() > 0);
    assert(Y.size() > 0);

    // Assuming regular grid
    unsigned int numSamples = X.at(0).size();

    for (unsigned int i = 0; i < numInputs; i++)
    {
        assert(X.at(i).size() == numSamples);
    }

    for (unsigned int i = 0; i < numOutputs; i++)
    {
        assert(Y.at(i).size() == numSamples);
    }

    /* Setup and solve equations Lc = R,
     * L = B'*W*B + l*D'*D
     * R = B'*W*y
     * c = control coefficients or knot averages.
     * B = basis functions at sample x-values,
     * W = weighting matrix for interpolating specific points
     * D = second-order finite difference matrix
     * l = penalizing parameter (increase for more smoothing)
     * y = sample y-values when calculating control coefficients,
     * y = sample x-values when calculating knot averages
     */

    SparseMatrix L, B;
    DenseMatrix Rx, Ry;

    SparseMatrix W;
    W.resize(numSamples, numSamples);
    W.setIdentity();

    controlPointEquationPsplineLHS(X, L, B, W, 0.03); //0.03
    controlPointEquationPsplineRHS(Y, Ry, B, W);
    controlPointEquationPsplineRHS(X, Rx, B, W);

    DenseMatrix CY;
    DenseMatrix CX;

    int numEquations = L.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if (!solveAsDense)
    {
        cout << "Computing B-spline control points using sparse solver." << endl;
        bool successfulSolve = (solveSparseLSE(L,Rx,CX) && solveSparseLSE(L,Ry,CY));

        solveAsDense = !successfulSolve;
    }

    if (solveAsDense)
    {
        cout << "Computing B-spline control points using dense solver." << endl;
        DenseMatrix Ld = L.toDense();
        bool successfulSolve = (solveDenseLSE(Ld,Rx,CX) && solveDenseLSE(Ld,Ry,CY));
        assert(successfulSolve);
    }

    coefficients = CY.transpose();
    knotaverages = CX.transpose();
}

bool Bspline::solveSparseLSE(const SparseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
{
    double errorTolerance = 1e-12;

    if (A.rows() != b.rows())
    {
        cout << "solveSparseLSE: incompatible matrix sizes" << endl;
        return false;
    }

    // Init BiCGSTAB solver
//    Eigen::BiCGSTAB<SparseMatrix> sparseSolver(A);

    // Init SparseLU solver
    Eigen::SparseLU<SparseMatrix > sparseSolver;
    // Compute the ordering permutation vector from the structural pattern of A
    sparseSolver.analyzePattern(A);
    // Compute the numerical factorization
    sparseSolver.factorize(A);

    if(sparseSolver.info() == Eigen::Success)
    {
        // Solve LSE
        x = sparseSolver.solve(b);

        if(!(sparseSolver.info() == Eigen::Success))
        {
            cout << "Sparse solver did not converge." << endl;
            return false;
        }

        double relativeError = (A*x - b).norm() / b.norm();
        if(relativeError > errorTolerance)
        {
            cout << "Sparse solver did not converge to acceptable tolerance." << endl;
            return false;
        }

        return true;
    }

    return false;
}

bool Bspline::solveDenseLSE(const DenseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const
{
    double errorTolerance = 1e-12;

    if (A.rows() != b.rows())
    {
        cout << "solveDenseLSE: incompatible matrix sizes" << endl;
        return false;
    }

    // Solve LSE
    x = A.colPivHouseholderQr().solve(b);

    double relativeError = (A*x - b).norm() / b.norm();
    if (relativeError > errorTolerance)
    {
        cout << "Dense solver did not converge to acceptable tolerance." << endl;
        return false;
    }

    return true;
}

void Bspline::controlPointEquationLHS(std::vector< std::vector<double> > &X, SparseMatrix &A)
{
    int vars = X.size();
    int samples = X.at(0).size();
    TensorIndex ti(basis.getTensorIndexDimension());

    int nnzPrCol = basis.supportedPrInterval();

    A.resize(ti.vectorSize(), samples); // Should be samples x vectorSize()
    A.reserve(DenseVector::Constant(samples, nnzPrCol));

    for (int i = 0; i < samples; i++)
    {
        DenseVector xi(vars);
        for (int j = 0; j < vars; j++)
        {
            xi(j) = X.at(j).at(i);
        }

        SparseVector tensorvaluesS = basis.evaluate(xi);

        for (SparseVector::InnerIterator it(tensorvaluesS); it; ++it)
        {
            A.insert(i,it.index()) = it.value();
        }
    }

    A.makeCompressed();
}

void Bspline::controlPointEquationRHS(std::vector< std::vector<double> > &Y, DenseMatrix &B)
{
    unsigned int vars = Y.size();
    unsigned int samples = Y.at(0).size();

    B.resize(samples, vars);

    for (unsigned int i = 0; i < vars; i++)
    {
        for (unsigned int j = 0; j < samples; j++)
        {
            B(j,i) = Y.at(i).at(j);
        }
    }
}

void Bspline::controlPointEquationPsplineLHS(std::vector< std::vector<double> > &X, SparseMatrix &L, SparseMatrix &B, SparseMatrix &W, double lambda)
{
    SparseMatrix D;

    // Basis function matrix
    getBasisFunctionMatrix(basis, X, B);

    // Second order finite difference matrix
    getSecondOrderFiniteDifferenceMatrix(D);

    L = B.transpose()*W*B + lambda*D.transpose()*D;
}

void Bspline::controlPointEquationPsplineRHS(std::vector< std::vector<double> > &Y, DenseMatrix &R, SparseMatrix &B, SparseMatrix &W)
{
    unsigned int vars = Y.size();
    unsigned int samples = Y.at(0).size();

    DenseMatrix C;
    C.resize(samples, vars);

    for (unsigned int i = 0; i < vars; i++)
    {
        for (unsigned int j = 0; j < samples; j++)
        {
            C(j,i) = Y.at(i).at(j);
        }
    }

    R = B.transpose()*W*C;
}

void Bspline::getBasisFunctionMatrix(Basis b, std::vector< std::vector<double> > &X, SparseMatrix &B)
{
    int vars = X.size();
    int samples = X.at(0).size();
    TensorIndex ti(b.getTensorIndexDimension());

    int nnzPrCol = b.supportedPrInterval();
    //B.resize(ti.vectorSize(), samples);
    B.resize(samples, ti.vectorSize());
    B.reserve(DenseVector::Constant(ti.vectorSize(), nnzPrCol)); // Want to reserve nnz per row, not col

    for (int i = 0; i < samples; i++)
    {
        DenseVector xi(vars);
        for (int j = 0; j < vars; j++)
        {
            xi(j) = X.at(j).at(i);
        }

        SparseVector tensorvaluesS = b.evaluate(xi);

        for (SparseVector::InnerIterator it(tensorvaluesS); it; ++it)
        {
            B.insert(i,it.index()) = it.value();
        }
    }

    B.makeCompressed();
}

// Function for generating second order finite-difference matrix, which is used for penalizing the
// (approximate) second derivative in control point calculation for P-splines.
void Bspline::getSecondOrderFiniteDifferenceMatrix(SparseMatrix &D)
{

    // Number of (total) basis functions - defines the number of columns in D
    int numCols = basis.numBasisFunctions();

    // Number of basis functions (and coefficients) in each dimension
    std::vector < int > dims = basis.getTensorIndexDimension();
    std::reverse(dims.begin(), dims.end());

    // Number of variables
    int vars = dims.size();

    for (int i=0; i < vars; i++)
    {
        // Need at least three coefficients in each dimension
        assert(basis.numBasisFunctions(i) >= 3);
    }

    // Number of rows in D and in each block
    int numRows = 0;
    std::vector< int > numBlkRows;
    for (int i = 0; i < vars; i++)
    {
        int prod = 1;
        for (int j = 0; j < vars; j++)
        {
            if (i == j)
                prod *= (dims[j] - 2);
            else
                prod *= dims[j];
        }
        numRows += prod;
        numBlkRows.push_back(prod);
    }

    // Resize and initialize D
    D.resize(numRows, numCols);                         // Resize (transpose because of reservation fn)
    D.reserve(DenseVector::Constant(numCols,2*vars));   // D has no more than two elems per col per dim

    int i = 0;                                          // Row index
    // Loop though each dimension (each dimension has its own block)
    for (int d = 0; d < vars; d++)
    {
        // Calculate left and right products
        int leftProd = 1;
        int rightProd = 1;
        for (int k = 0; k < d; k++)
        {
            leftProd *= dims[k];
        }
        for (int k = d+1; k < vars; k++)
        {
            rightProd *= dims[k];
        }

        // Loop through subblocks on the block diagonal
        for (int j = 0; j < rightProd; j++)
        {
            // Start column of current subblock
            int blkBaseCol = j*leftProd*dims[d];
            // Block rows [I -2I I] of subblock
            for (int l = 0; l < (dims[d] - 2); l++)
            {
                // Special case for first dimension
                if (d == 0)
                {
                    int k = j*leftProd*dims[d] + l;
                    D.insert(i,k) = 1;
                    k += leftProd;
                    D.insert(i,k) = -2;
                    k += leftProd;
                    D.insert(i,k) = 1;
                    i++;
                }
                else
                {
                    // Loop for identity matrix
                    for (int n = 0; n < leftProd; n++)
                    {
                        int k = blkBaseCol + l*leftProd + n;
                        D.insert(i,k) = 1;
                        k += leftProd;
                        D.insert(i,k) = -2;
                        k += leftProd;
                        D.insert(i,k) = 1;
                        i++;
                    }
                }
            }
        }
    }

//    // Number of (total) basis functions - defines the number of columns in D
//    int numCols = basis.numBasisFunctions();

//    // Number of basis functions (and coefficients) in each dimension
//    std::vector < int > dims = basis.getTensorIndexDimension();
//    std::reverse(dims.begin(), dims.end()); // flip vector

//    // Number of variables
//    int vars = dims.size();

//    for (int i=0; i < vars; i++)
//    {
//        // Need at least three coefficients in each dimension
//        assert(basis.numBasisFunctions(i) >= 3);
//        dims.push_back(basis.numBasisFunctions(i));
//    }

//    // Number of rows in D and in each block
//    int numRows = 0;
//    std::vector< int > numBlkRows;
//    for (int i = 0; i < vars; i++)
//    {
//        int prod = 1;
//        for (int j = 0; j < vars; j++)
//        {
//            if (i == j)
//                prod *= (dims[j] - 2);
//            else
//                prod *= dims[j];
//        }
//        numRows += prod;
//        numBlkRows.push_back(prod);
//    }

//    // Resize and initialize D
//    DenseMatrix Dd;
//    Dd.resize(numRows, numCols);
//    Dd.setZero();
//    cout << numRows << "/" << numCols << endl;

//    // Accumulate variable for product of preceding dimensions
//    int acc = 1;

//    // Accumulate variable for block insertion point
//    int insertionRow = 0;

//    // Base block
//    DenseMatrix d(1,3);
//    d << 1, -2, 1;

//    // Calculate block matrix associated with each variable and insert into D
//    for (int i = 0; i < vars; i++)
//    {

//        DenseMatrix blk;
//        DenseMatrix subBlk(acc*(dims[i]-2),acc*dims[i]);  subBlk.setZero();
//        for (int j = 0; j < dims[i]-2; j++)
//        {
//            DenseMatrix Ij(acc, acc); Ij.setIdentity();
//            subBlk.block(acc*j,acc*j,acc,3*acc) = kroneckerProduct(d, Ij);
//        }
//        // Complete Kronecker products
//        DenseMatrix tmp1, tmp2;
//        tmp2 = subBlk;

//        for (int j = i+1; j < vars; j++)
//        {
//            tmp1 = tmp2;
//            DenseMatrix Ij(dims[j],dims[j]); Ij.setIdentity();
//            tmp2 = kroneckerProduct(Ij, tmp1);
//        }
//        blk = tmp2;
//        acc *= dims[i];

//        // Insert block into D matrix
//        if (i > 0) insertionRow += numBlkRows[i-1];
//        Dd.block(insertionRow,0,numBlkRows[i],numCols) = blk;
//    }

//    // Convert to sparse matrix and compress
//    D = Dd.sparseView();

    D.makeCompressed();
}

bool Bspline::insertKnots(double tau, unsigned int dim, unsigned int multiplicity)
{
    // Test multiplicity at knot
    if (basis.getKnotMultiplicity(dim, tau) + multiplicity > basis.getBasisDegree(dim) + 1)
        return false;

    // Insert knots and compute knot insertion matrix
    SparseMatrix A;
    if (!basis.insertKnots(A, tau, dim, multiplicity))
        return false;

    // Update control points
//    DenseMatrix P = getControlPoints(); // TODO: remove unnecessary copy
//    assert(A.cols() == P.cols());
//    DenseMatrix newP = P*A.transpose();
//    setControlPoints(newP);

    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();

    return true;
}

bool Bspline::refineKnotSequences()
{
    // Compute knot insertion matrix
    SparseMatrix A;
    if (!basis.refineKnots(A))
        return false;

    // Update control points
//    DenseMatrix P = getControlPoints(); // TODO: remove unnecessary copy
//    assert(A.cols() == P.cols());
//    DenseMatrix newP = P*A.transpose();
//    setControlPoints(newP);

    assert(A.cols() == coefficients.cols());
    coefficients = coefficients*A.transpose();
    knotaverages = knotaverages*A.transpose();

    return true;
}

// NOTE: Do this lower in the hierarchy so that all knots can be added at once
bool Bspline::regularSequences(std::vector<double> &lb, std::vector<double> &ub)
{
//    assert(checkBsplinedata(knotsequences, controlpoints));
    // Add and remove controlpoints and knots to make the b-spline p-regular with support [lb, ub]
    if (!(lb.size() == numInputs && ub.size() == numInputs))
        return false;

    for (unsigned int dim = 0; dim < numInputs; dim++)
    {
        unsigned int multiplicityTarget = basis.getBasisDegree(dim) + 1;

        // Inserting many knots at the time (to save number of B-spline coefficient calculations)
        // NOTE: This method generates knot insertion matrices with more nonzero elements than
        // the method that inserts one knot at the time. This causes the preallocation of
        // kronecker product matrices to become too small and the speed deteriorates drastically
        // in higher dimensions because reallocation is necessary. This can be prevented by
        // precomputing the number of nonzeros when preallocating memory (see myKroneckerProduct).
        int numKnotsLB = multiplicityTarget - basis.getKnotMultiplicity(dim, lb.at(dim));
        if (numKnotsLB > 0)
        {
            if (!insertKnots(lb.at(dim), dim, numKnotsLB))
                return false;
        }

        int numKnotsUB = multiplicityTarget - basis.getKnotMultiplicity(dim, ub.at(dim));
        if (numKnotsUB > 0)
        {
            if (!insertKnots(ub.at(dim), dim, numKnotsUB))
                return false;
        }

        // Old insertion method: inserts one knot at the time
//        while (basis.getKnotMultiplicity(dim, lb.at(dim)) < multiplicityTarget)
//        {
//            insertKnots(lb.at(dim), dim);
//        }

//        while (basis.getKnotMultiplicity(dim, ub.at(dim)) < multiplicityTarget)
//        {
//            insertKnots(ub.at(dim), dim);
//        }
    }

//    assert(checkBsplinedata(knotsequences, controlpoints));
    return true;
}

bool Bspline::removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub)
{
    assert(lb.size() == numInputs);
    assert(ub.size() == numInputs);

    SparseMatrix A;
    if (!basis.reduceSupport(lb, ub, A))
        return false;

    if (coefficients.cols() != A.rows())
        return false;

    // Remove unsupported control points (basis functions)
    coefficients = coefficients*A;
    knotaverages = knotaverages*A;

    return true;
}

//// Tests BSpline at 100 randomly chosen points
//bool BSpline::testBspline(DenseVector (*afunc)(const DenseVector))
//{
//    // Sample a little bit
//    double reltol = 1e-6;

//    DenseVector x; x.setZero(numInputs);

//    for (int i = 0; i < 100; i++)
//    {
//        // Generate a random number
//        for (unsigned int j = 0; j < numInputs; j++)
//        {
//            x(j) = randomInteger(domainLowerBound.at(j),domainUpperBound.at(j));
//        }

//        // Evaluate function
//        DenseVector y = afunc(x);
//        assert(y.rows() == numOutputs);

//        // Evaluate bspline
//        DenseVector ys = this->evaluate(x);

//        // Evaluate error
//        DenseVector e = y-ys;
//        if (y.norm() != 0) e = e/y.norm();
//        if (e.norm() > reltol) return false;
//    }

//    return true;
//}

#ifndef TENSORPRODUCTBSPLINE_H
#define TENSORPRODUCTBSPLINE_H

#include "generaldefinitions.h"
//#include "../interpolationtable.h"
#include "sorteddatatable.h"
#include "tensorindex.h"
#include "basis.h"

// Enum for different B-spline types
enum class BsplineType
{
    EXPLICIT,   // B-spline is explicitly given by control points and knot sequences
    FREE,       // Interpolates all points. Ensures p-1'th derivative continuous at x(2) and x(n-1). p+1-regular knot sequence with two deleted knots
    NATURAL,    // Not implemented. Interpolates all points. Ensures second derivative of B-spline is zero at end points.
    PSPLINE     // Minimizes objective which penalizes both deviation (for interpolation) and second derivative (for smoothing)
};

// Maybe change name to EXPLICIT, CUBIC_FREE, CUBIC_NATURAL, PSPLINE,

// Implements tensor product B-splines
class Bspline
{
public:

    // Construct B-spline from knot sequences, control coefficients (assumed vectorized), and basis degrees
    //Bspline(std::vector<double> coefficients, std::vector<double> knotSequence, int basisDegrees);
    //Bspline(std::vector<double> coefficients, std::vector< std::vector<double> > knotSequences, std::vector<int> basisDegrees);
    Bspline(DenseMatrix coefficients, std::vector< std::vector<double> > knotSequences, std::vector<int> basisDegrees);

    // Construct B-spline from interpolation data (InterpolationTable)
//    Bspline(InterpolationTable &data, int basisdegree);
//    Bspline(InterpolationTable &data, std::vector<int> basisdegree);

    // Construct B-spline from interpolation data (SortedDataTable)
    Bspline(SortedDataTable &data, int basisdegree);
    Bspline(SortedDataTable &data, std::vector<int> basisdegree);

    virtual Bspline* clone() const { return new Bspline(*this); }

    void init();

    // Evaluation of B-spline
    DenseVector evaluate(DenseVector &x);
    DenseMatrix jacobian(DenseVector &x);
    DenseMatrix hessian(DenseVector &x); // Supports only 1 output

    // Getters
    unsigned int getNumInputs() const { return numInputs; }
    unsigned int getNumOutputs() const { return numOutputs; }
    unsigned int getNumControlPoints() const { return coefficients.cols(); }

    std::vector< std::vector<double> > getKnotVectors() const;

    std::vector<double> getDomainUpperBound() const;
    std::vector<double> getDomainLowerBound() const;

    // Control point related
    void setControlPoints(DenseMatrix &controlPoints);
    DenseMatrix getControlPoints() const;
    bool checkControlPoints() const;

    // B-spline operations
    bool reduceDomain(std::vector<double> lb, std::vector<double> ub, bool regularKnotsequences = true, bool refineKnotsequences = true);

    bool insertKnots(double tau, unsigned int dim, unsigned int multiplicity = 1); // MOVE BACK TO PRIVATE

private:

    Basis basis;
    DenseMatrix knotaverages; // One row per input
    DenseMatrix coefficients; // One row per output

    unsigned int numInputs; // Dimension of x
    unsigned int numOutputs; // Dimension of y = f(x)

    // Control point calculations
    void calculateKnotAverages();
    void calculateControlPoints(std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y);
    void controlPointEquationLHS(std::vector< std::vector<double> > &X, SparseMatrix &A);
    void controlPointEquationRHS(std::vector< std::vector<double> > &Y, DenseMatrix &B);
    bool solveSparseLSE(const SparseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const;
    bool solveDenseLSE(const DenseMatrix &A, const DenseMatrix &b, DenseMatrix &x) const;

    // P-spline control point calculation
    void calculateControlPointsPspline(std::vector< std::vector<double> > &X, std::vector< std::vector<double> > &Y);
    void controlPointEquationPsplineLHS(std::vector< std::vector<double> > &X, SparseMatrix &L, SparseMatrix &B, SparseMatrix &W, double lambda);
    void controlPointEquationPsplineRHS(std::vector< std::vector<double> > &Y, DenseMatrix &R, SparseMatrix &B, SparseMatrix &W);
    void getSecondOrderFiniteDifferenceMatrix(SparseMatrix &D);
    void getBasisFunctionMatrix(Basis b, std::vector< std::vector<double> > &X, SparseMatrix &B);

    // Domain reduction
    bool regularSequences(std::vector<double> &lb, std::vector<double> &ub);
    bool removeUnsupportedBasisFunctions(std::vector<double> &lb, std::vector<double> &ub);

    // Knot insertion and refinement
    bool refineKnotSequences(); // All knots in one shabang

    // Helper functions
    bool valueInsideDomain(DenseVector x);

    friend class TBtestbench;
};

#endif // TENSORPRODUCTBSPLINE_H

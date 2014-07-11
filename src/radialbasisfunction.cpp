#include "radialbasisfunction.h"
#include "linearsolvers.h"

using std::cout;
using std::endl;

RadialBasisFunction::RadialBasisFunction(SortedDataTable &samples, RadialBasisFunctionType type)
    : RadialBasisFunction(samples, type, true)
{
}

RadialBasisFunction::RadialBasisFunction(SortedDataTable &samples, RadialBasisFunctionType type, bool normalized)
    : samples(samples), normalized(normalized), dim(samples.getDimX()), numSamples(samples.getNumSamples())
{
    if (type == RadialBasisFunctionType::THIN_PLATE_SPLINE)
    {
        fn = std::shared_ptr<RBF_fn>(new ThinPlateSpline());
    }
    else if (type == RadialBasisFunctionType::MULTIQUADRIC)
    {
        fn = std::shared_ptr<RBF_fn>(new Multiquadric());
    }
    else if (type == RadialBasisFunctionType::INVERSE_QUADRIC)
    {
        fn = std::shared_ptr<RBF_fn>(new InverseQuadric());
    }
    else if (type == RadialBasisFunctionType::INVERSE_MULTIQUADRIC)
    {
        fn = std::shared_ptr<RBF_fn>(new InverseMultiquadric());
    }
    else if (type == RadialBasisFunctionType::GAUSSIAN)
    {
        fn = std::shared_ptr<RBF_fn>(new Gaussian());
    }
    else
    {
        fn = std::shared_ptr<RBF_fn>(new ThinPlateSpline());
    }

    /* Want to solve the linear system A*w = b,
     * where w is the vector of weights.
     */
    //SparseMatrix A(numSamples,numSamples);
    //A.reserve(numSamples*numSamples);
    DenseMatrix A; A.setZero(numSamples,numSamples);
    DenseMatrix b; b.setZero(numSamples,1);
    cout << "Starting loop" << endl;
    int i=0;
    std::multiset<DataSample>::const_iterator it1, it2;
    for (it1 = samples.cbegin(); it1 != samples.cend(); ++it1, ++i)
    {
        double sum = 0;
        int j=0;
        for (it2 = samples.cbegin(); it2 != samples.cend(); ++it2, ++j)
        {
            double val = fn->eval(dist((*it1).getX(), (*it2).getX()));
            if (val != 0)
            {
                //A.insert(i,j) = val;
                A(i,j) = val;
                sum += val;
            }
        }

        double val = (*it1).getY().at(0);
        if (normalized) b(i) = sum*val;
        else b(i) = val;
    }

    //A.makeCompressed();

    int numEquations = A.rows();
    int maxNumEquations = pow(2,10);

    bool solveAsDense = (numEquations < maxNumEquations);

    if (!solveAsDense)
    {
        cout << "Computing RBF weights using sparse solver." << endl;
        SparseMatrix As = A.sparseView();
        SparseLU s;
        solveAsDense = !s.solve(As,b,weights);
    }

    if (solveAsDense)
    {
        cout << "Computing RBF weights using dense solver." << endl;
        //DenseMatrix Ad = A.toDense();
        DenseQR s;
        bool success = s.solve(A,b,weights);
        assert(success);
    }
}

double RadialBasisFunction::evaluate(DenseVector &x) const
{
    std::vector<double> y;
    for (int i=0;i<x.rows();i++)
        y.push_back(x(i));
    return evaluate(y);
}

double RadialBasisFunction::evaluate(std::vector<double> &x) const
{
    assert(x.size() == dim);
    double fval, sum=0, sumw=0;
    int i = 0;
    std::multiset<DataSample>::const_iterator it;
    for (it = samples.cbegin(); it != samples.cend(); ++it, ++i)
    {
        fval = fn->eval(dist(x,(*it).getX()));
        sumw += weights(i)*fval;
        sum += fval;
    }
    return normalized ? sumw/sum : sumw;
}


/*
 * Computes Euclidean distance ||x-y||
 */
double RadialBasisFunction::dist(std::vector<double> x, std::vector<double> y) const
{
    assert(x.size() == y.size());
    double sum = 0;
    for (unsigned int i=0; i<x.size(); i++)
        sum += (x.at(i)-y.at(i))*(x.at(i)-y.at(i));
    return std::sqrt(sum);
}

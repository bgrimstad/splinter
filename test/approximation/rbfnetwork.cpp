/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <Catch.h>
#include <testingutilities.h>
#include "rbfnetwork.h"
#include <testfunction.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][rbf][radialbasisfunction]"
#define COMMON_TEXT " value approximation test with polynomials"


/*
 * Gaussian
 */
TEST_CASE("Gaussian RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[RBFType::gaussian][function-value]")
{
    double one_eps = 2.1e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::GAUSSIAN);
                             },
                             400,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Gaussian RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[RBFType::gaussian][jacobian]")
{
    double one_eps = 1.5e-1;
    double two_eps = 1.5e-1;
    double inf_eps = 1.5e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::GAUSSIAN);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Gaussian RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[RBFType::gaussian][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RBFType::GAUSSIAN);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Inverse multiquadric
 */
TEST_CASE("Inverse Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-multiquadric][function-value]")
{
    double one_eps = 2e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::INVERSE_MULTIQUADRIC);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Inverse Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-multiquadric][jacobian]")
{
    double one_eps = 1e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::INVERSE_MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Inverse Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-multiquadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RBFType::INVERSE_MULTIQUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Inverse quadric
 */
TEST_CASE("Inverse Quadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-quadric][function-value]")
{
    double one_eps = 2.2e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::INVERSE_QUADRIC);
                             },
                             350,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Inverse Quadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-quadric][jacobian]")
{
    double one_eps = 1e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::INVERSE_QUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Inverse Quadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[RBFType::inverse-quadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RBFType::INVERSE_QUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Multiquadric
 */
TEST_CASE("Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[RBFType::multiquadric][function-value]")
{
    double one_eps = 2.2e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::MULTIQUADRIC);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[RBFType::multiquadric][jacobian]")
{
    double one_eps = 1e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[RBFType::multiquadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RBFType::MULTIQUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Thin plate spline
 */
TEST_CASE("Thin plate spline RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[RBFType::thin-plate-spline][function-value]")
{
    double one_eps = 2e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::THIN_PLATE_SPLINE);
                             },
                             300,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Thin plate spline RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[RBFType::thin-plate-spline][jacobian]")
{
    double one_eps = 1e-1;
    double two_eps = 1e-1;
    double inf_eps = 1e-1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Function *) new RBFNetwork(table, RBFType::THIN_PLATE_SPLINE);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Thin plate spline RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[RBFType::thin-plate-spline][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RBFType::THIN_PLATE_SPLINE);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}

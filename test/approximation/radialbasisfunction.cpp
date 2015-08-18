/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "term.h"
#include <Catch.h>
#include <testingutilities.h>
#include <radialbasisfunction.h>

using namespace SPLINTER;


#define COMMON_TAGS "[approximation][polynomial][rbf][radialbasisfunction]"
#define COMMON_TEXT " value approximation test with polynomials"


/*
 * Gaussian
 */
TEST_CASE("Gaussian RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][function-value]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
                             },
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Gaussian RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][jacobian]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Gaussian RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::gaussian][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::GAUSSIAN);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Inverse multiquadric
 */
TEST_CASE("Inverse Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][function-value]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Inverse Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][jacobian]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Inverse Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-multiquadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_MULTIQUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Inverse quadric
 */
TEST_CASE("Inverse Quadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][function-value]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
                             },
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Inverse Quadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][jacobian]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Inverse Quadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::inverse-quadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::INVERSE_QUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Multiquadric
 */
TEST_CASE("Multiquadric RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][function-value]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Multiquadric RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][jacobian]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Multiquadric RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::multiquadric][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::MULTIQUADRIC);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}


/*
 * Thin plate spline
 */
TEST_CASE("Thin plate spline RadialBasisFunction function" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][function-value]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareFunctionValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
                             },
                             200,  // Number of points to sample at
                             1337, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

TEST_CASE("Thin plate spline RadialBasisFunction jacobian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][jacobian]")
{
    double one_eps = 0.1;
    double two_eps = 0.1;
    double inf_eps = 0.1;

    for(auto testFunc : getPolynomialFunctions())
    {
        compareJacobianValue(testFunc,
                             [](const DataTable &table)
                             {
                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
                             },
                             200,  // Number of points to sample at
                             999, // Number of points to test against
                             one_eps, two_eps, inf_eps);
    }
}

// TODO: Uncomment when implemented
//TEST_CASE("Thin plate spline RadialBasisFunction hessian" COMMON_TEXT, COMMON_TAGS "[radialbasisfunctiontype::thin-plate-spline][hessian]")
//{
//    for(auto testFunc : getPolynomialFunctions())
//    {
//        checkHessianSymmetry(testFunc,
//                             [](const DataTable &table)
//                             {
//                                 return (Approximant *) new RadialBasisFunction(table, RadialBasisFunctionType::THIN_PLATE_SPLINE);
//                             },
//                             200,   // Number of points to sample at
//                             999); // Number of points to test against
//    }
//}

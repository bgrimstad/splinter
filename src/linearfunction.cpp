/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "linearfunction.h"

namespace SPLINTER
{

//double LinearFunction::eval(DenseVector x) const
//{
//    if (x.size() != numVariables)
//        throw Exception("LinearFunction::eval: Wrong dimension on evaluation point x.");
//    auto res = coefficients.transpose()*evalBasis(x);
//    return res(0);
//}
//
//DenseMatrix LinearFunction::evalJacobian(DenseVector x) const
//{
//    if (x.size() != numVariables)
//        throw Exception("LinearFunction::eval: Wrong dimension on evaluation point x.");
//    return coefficients.transpose()*evalBasisJacobian(x);
//}

} // namespace SPLINTER

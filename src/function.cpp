/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <function.h>
#include "utilities.h"

namespace SPLINTER
{

DenseVector Function::eval(const DenseVector &x) const {
    return std_to_eig_vec(eval(eig_to_std_vec(x)));
}

std::vector<std::vector<double>> Function::eval_jacobian(const std::vector<double> &x) const
{
    auto denseX = std_to_eig_vec(x);

    return eig_mat_to_std_vec_vec(eval_jacobian(denseX));
}

DenseMatrix Function::eval_jacobian(const DenseVector &x) const
{
    return central_difference(x);
}

DenseMatrix Function::central_difference(const DenseVector &x) const
{
    DenseMatrix Jac(dim_y, dim_x);

    double h = 1e-6; // perturbation step size
    double hForward = 0.5*h;
    double hBackward = 0.5*h;

    for (unsigned int i = 0; i < dim_x; ++i)
    {
        DenseVector xForward(x);
        xForward(i) = xForward(i) + hForward;

        DenseVector xBackward(x);
        xBackward(i) = xBackward(i) - hBackward;

        auto yForward = eval(xForward);
        auto yBackward = eval(xBackward);

        for (unsigned int j = 0; j < dim_y; ++j)
            Jac(j, i) = (yForward(j) - yBackward(j)) / (hBackward + hForward);
    }

    return Jac;
}

void Function::check_input(const DenseVector &x) const {
    return check_input(eig_to_std_vec(x));
}


} // namespace SPLINTER
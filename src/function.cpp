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
    auto eig_x = std_to_eig_vec(x);

    return eig_to_std_mat(eval_jacobian(eig_x));
}

DenseMatrix Function::eval_jacobian(const DenseVector &x) const
{
    return central_difference(x);
}

DenseMatrix Function::central_difference(const DenseVector &x) const
{
    DenseMatrix jacobian(dim_y, dim_x);

    double h = 1e-6; // perturbation step size
    double h_forward = 0.5*h;
    double h_backward = 0.5*h;

    for (unsigned int i = 0; i < dim_x; ++i)
    {
        DenseVector x_forward(x);
        x_forward(i) = x_forward(i) + h_forward;

        DenseVector x_backward(x);
        x_backward(i) = x_backward(i) - h_backward;

        auto y_forward = eval(x_forward);
        auto y_backward = eval(x_backward);

        for (unsigned int j = 0; j < dim_y; ++j)
            jacobian(j, i) = (y_forward(j) - y_backward(j)) / (h_backward + h_forward);
    }

    return jacobian;
}

void Function::check_input(const DenseVector &x) const {
    return check_input(eig_to_std_vec(x));
}


} // namespace SPLINTER
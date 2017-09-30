/*
 * This file is part of the SPLINTER library.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SPLINTER_KRONECKER_PRODUCT_H
#define SPLINTER_KRONECKER_PRODUCT_H

#include "definitions.h"

namespace SPLINTER
{

SparseMatrix my_kronecker_product(const SparseMatrix &A, const SparseMatrix &B);

// Apply Kronecker product on several vectors or matrices
SparseVector kronecker_product_vectors(const std::vector<SparseVector> &vectors);
DenseVector kronecker_product_vectors(const std::vector<DenseVector> &vectors);
SparseMatrix kronecker_product_matrices(const std::vector<SparseMatrix> &matrices);

} // namespace SPLINTER

#endif // SPLINTER_KRONECKER_PRODUCT_H

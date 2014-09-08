#ifndef MYKRONECKERPRODUCT_H
#define MYKRONECKERPRODUCT_H

#include "generaldefinitions.h"

namespace MultivariateSplines
{

void myKroneckerProduct(const SparseMatrix &A, const SparseMatrix &B, SparseMatrix &AB);

} // namespace MultivariateSplines

#endif // MYKRONECKERPRODUCT_H

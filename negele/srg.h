// Copyright 2020 Matthias Heinz
#ifndef NEGELE_SRG_H_
#define NEGELE_SRG_H_

#include "negele/matrix.h"

namespace negele {
namespace srg {
void rhs(
    const negele::matrix::MomentumSpaceMatrix& v,
    negele::matrix::MomentumSpaceMatrix& dvdt,  // NOLINT [runtime/references]
    const double lambda);
}
}  // namespace negele

#endif  // NEGELE_SRG_H_

// Copyright 2020 Matthias Heinz
#include "negele/srg.h"

#include <math.h>

#include <iostream>

#include "negele/matrix.h"

void negele::srg::rhs(
    const negele::matrix::MomentumSpaceMatrix& v,
    negele::matrix::MomentumSpaceMatrix& dvdt,  // NOLINT [runtime/references]
    const double lambda) {
  dvdt.clear_mels();
  auto trel = v.kinetic_energy();
  auto v_mels = v.matrix_elements_with_weights();
  auto ham = trel + v_mels;
  auto evs = ham.eigenvalues();
  std::cout << lambda << ", " << evs[0] << ", " << evs[1] << ", " << ham[30]
            << std::endl;
  auto dvdt_mels =
      -4.0 / pow(lambda, 5) *
      negele::matrix::commutator(negele::matrix::commutator(trel, ham), ham);
  dvdt.set_matrix_elements(dvdt_mels);
}

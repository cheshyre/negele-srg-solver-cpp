// Copyright 2020 Matthias Heinz
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "negele/matrix.h"
#include "negele/potential.h"
#include "negele/quad.h"
#include "negele/srg.h"

int main(void) {
  auto quadrature = negele::quad::get_gauss_legendre_quadrature(100, 0.0, 25.0);
  auto pts = std::get<0>(quadrature);
  auto weights = std::get<1>(quadrature);

  int dim = pts.size();

  double v1 = 12.0;
  double sig1 = 0.2;
  double v2 = -12.0;
  double sig2 = 0.8;

  negele::potential::NegelePotential pot(v1, v2, sig1, sig2);

  std::vector<double> pot_mels(dim * dim, 0.0);
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      pot_mels[i * dim + j] =
          pot.eval(pts[j], pts[i]) + pot.eval(pts[j], -1 * pts[i]);
    }
  }

  // for (int i = 0; i < dim; i++) {
  //   std::cout << pot_mels[i] << ", ";
  // }
  // std::cout << std::endl;

  // negele::matrix::SquareMatrix pot_matrix(dim, pot_mels);

  negele::matrix::MomentumSpaceMatrix pot_matrix(pts, weights, pot_mels);

  auto ham_matrix =
      pot_matrix.matrix_elements_with_weights() + pot_matrix.kinetic_energy();

  auto evs = ham_matrix.eigenvalues();
  for (int i = 0; i < 10; i++) {
    std::cout << evs[i] << std::endl;
  }

  // typedef boost::numeric::odeint::runge_kutta4<
  //     negele::matrix::MomentumSpaceMatrix, double,
  //     negele::matrix::MomentumSpaceMatrix, double,
  //     boost::numeric::odeint::vector_space_algebra>
  //     stepper_type;
  // stepper_type stepper;

  // const double dlambda = -0.01;
  // for (double lambda = 50.0; lambda >= 10.0; lambda += dlambda) {
  //   stepper.do_step(negele::srg::rhs, pot_matrix, lambda, dlambda);
  //   auto ham_matrix_new =
  //       pot_matrix.matrix_elements_with_weights() +
  //       pot_matrix.kinetic_energy();

  //   auto evs_new = ham_matrix_new.eigenvalues();
  //   std::cout << lambda << ", " << evs_new[0] << ", " << ham_matrix_new[30]
  //             << std::endl;
  // }

  typedef boost::numeric::odeint::runge_kutta_dopri5<
      negele::matrix::MomentumSpaceMatrix, double,
      negele::matrix::MomentumSpaceMatrix, double,
      boost::numeric::odeint::vector_space_algebra>
      stepper;
  int steps = boost::numeric::odeint::integrate_adaptive(
      boost::numeric::odeint::make_controlled<stepper>(1E-6, 1E-6),
      negele::srg::rhs, pot_matrix, 50.0, 2.0, -0.1);

  return 0;
}

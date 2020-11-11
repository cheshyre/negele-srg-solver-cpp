// Copyright 2020 Matthias Heinz
#ifndef NEGELE_QUAD_H_
#define NEGELE_QUAD_H_

#include <tuple>
#include <vector>

namespace negele {
namespace quad {

/**
 * @brief Get the gauss legendre quadrature object
 *
 * @param num_pts
 * @param x_min
 * @param x_max
 * @return std::tuple<std::vector, std::vector>
 */
std::tuple<std::vector<double>, std::vector<double>>
get_gauss_legendre_quadrature(int num_pts, double x_min, double x_max);

}  // namespace quad
}  // namespace negele

#endif  // NEGELE_QUAD_H_

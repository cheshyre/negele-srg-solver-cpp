// Copyright 2020 Matthias Heinz
#include "negele/matrix.h"

#include <lapack.h>

#include <cmath>
#include <iostream>
#include <vector>

#include <boost/numeric/odeint.hpp>

#include "negele/constants.h"

namespace boost {
namespace numeric {
namespace odeint {

template <>
struct is_resizeable<negele::matrix::MomentumSpaceMatrix> {  // declare
                                                             // resizeability
  typedef boost::true_type type;
  const static bool value = type::value;
};

template <>
struct same_size_impl<negele::matrix::MomentumSpaceMatrix,
                      negele::matrix::MomentumSpaceMatrix> {  // define how to
                                                              // check size
  static bool same_size(const negele::matrix::MomentumSpaceMatrix& v1,
                        const negele::matrix::MomentumSpaceMatrix& v2) {
    return v1.same_size(v2);
  }
};

template <>
struct resize_impl<negele::matrix::MomentumSpaceMatrix,
                   negele::matrix::MomentumSpaceMatrix> {  // define how to
                                                           // resize
  static void resize(negele::matrix::MomentumSpaceMatrix& v1,
                     const negele::matrix::MomentumSpaceMatrix& v2) {
    v1.resize(v2);
  }
};

}  // namespace odeint
}  // namespace numeric
}  // namespace boost

negele::matrix::SquareMatrix::SquareMatrix(int dim, double value)
    : m_dim(dim), m_matrix(dim * dim, value) {}

negele::matrix::SquareMatrix::SquareMatrix(int dim,
                                           const std::vector<double>& mels)
    : m_dim(dim), m_matrix(mels) {
  if (m_matrix.size() != m_dim * m_dim) {
    std::cout << "Matrix initialized with incompatible dimensions: \n"
              << "\tsingle-axis dimension given: " << m_dim << "\n"
              << "\tactual size of matrix: " << m_matrix.size() << "\n";
    exit(EXIT_FAILURE);
  }
}

negele::matrix::SquareMatrix negele::matrix::SquareMatrix::from_diag(
    const std::vector<double>& diag_mels) {
  int dim = diag_mels.size();
  std::vector<double> mels(dim * dim, 0.0);
  for (int i = 0; i < dim; i++) {
    mels[i * dim + i] = diag_mels[i];
  }

  return negele::matrix::SquareMatrix(dim, mels);
}

negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::operator+=(
    const negele::matrix::SquareMatrix& other) {
  if (m_dim != other.m_dim) {
    std::cout << "Matrices have incompatible dimensions: " << m_dim << ", and "
              << other.m_dim << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_matrix[i] += other.m_matrix[i];
  }
  return *this;
}

negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::operator-=(
    const negele::matrix::SquareMatrix& other) {
  if (m_dim != other.m_dim) {
    std::cout << "Matrices have incompatible dimensions: " << m_dim << ", and "
              << other.m_dim << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_matrix[i] -= other.m_matrix[i];
  }
  return *this;
}
negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::operator*=(
    double factor) {
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_matrix[i] *= factor;
  }
  return *this;
}
negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::operator/=(
    double factor) {
  if (factor == 0.0) {
    std::cout << "Attempted division by 0.\n";
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_matrix[i] /= factor;
  }
  return *this;
}

negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::operator/=(
    const negele::matrix::SquareMatrix& other) {
  if (m_dim != other.m_dim) {
    std::cout << "Matrices have incompatible dimensions: " << m_dim << ", and "
              << other.m_dim << "\n";
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < m_dim * m_dim; i++) {
    if (other.m_matrix[i] != 0.0) {
      m_matrix[i] /= other.m_matrix[i];
    }
  }
  return *this;
}

negele::matrix::SquareMatrix& negele::matrix::SquareMatrix::abs() {
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_matrix[i] = fabs(m_matrix[i]);
  }
  return *this;
}

double negele::matrix::SquareMatrix::max() const {
  double val = 0.0;
  if (m_dim > 0) {
    val = m_matrix[0];
  }
  for (int i = 0; i < m_dim * m_dim; i++) {
    if (m_matrix[i] > val) {
      val = m_matrix[i];
    }
  }
  return val;
}

negele::matrix::SquareMatrix negele::matrix::operator+(
    const SquareMatrix& a, const negele::matrix::SquareMatrix& b) {
  negele::matrix::SquareMatrix c(a);
  c += b;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::operator-(
    const negele::matrix::SquareMatrix& a,
    const negele::matrix::SquareMatrix& b) {
  negele::matrix::SquareMatrix c(a);
  c -= b;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::operator*(
    double factor, const negele::matrix::SquareMatrix& mat) {
  negele::matrix::SquareMatrix c(mat);
  c *= factor;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::operator*(
    const negele::matrix::SquareMatrix& mat, double factor) {
  negele::matrix::SquareMatrix c(mat);
  c *= factor;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::operator/(
    const negele::matrix::SquareMatrix& mat, double factor) {
  negele::matrix::SquareMatrix c(mat);
  c /= factor;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::operator/(
    const negele::matrix::SquareMatrix& a,
    const negele::matrix::SquareMatrix& b) {
  negele::matrix::SquareMatrix c(a);
  c /= b;
  return c;
}
negele::matrix::SquareMatrix negele::matrix::abs(
    const negele::matrix::SquareMatrix& a) {
  negele::matrix::SquareMatrix c(a);
  c.abs();
  return c;
}

negele::matrix::SquareMatrix negele::matrix::operator*(
    const negele::matrix::SquareMatrix& a,
    const negele::matrix::SquareMatrix& b) {
  if (a.m_dim != b.m_dim) {
    std::cout << "Matrices have incompatible dimensions: " << a.m_dim
              << ", and " << b.m_dim << "\n";
    exit(EXIT_FAILURE);
  }
  auto new_mels = negele::matrix::mat_mul(a.m_matrix, b.m_matrix, a.m_dim);
  return negele::matrix::SquareMatrix(a.m_dim, new_mels);
}

std::vector<double> negele::matrix::SquareMatrix::eigenvalues() const {
  std::vector<double> evs(m_dim, 0.0);
  std::vector<double> matrix_copy(m_matrix);
  std::vector<double> work(20 * m_dim, 0.0);
  int info = 0;
  char uplo = 'U';
  char jobz = 'N';
  int lwork = work.size();
  LAPACK_dsyev(&jobz, &uplo, &m_dim, matrix_copy.data(), &m_dim, evs.data(),
               work.data(), &lwork, &info);

  return evs;
}

negele::matrix::MomentumSpaceMatrix::MomentumSpaceMatrix(
    const std::vector<double>& pts, const std::vector<double>& weights,
    const std::vector<double>& mels)
    : m_dim(pts.size()),
      m_pts(pts),
      m_weights(weights),
      m_mels_with_weights(pts.size(), 0.0),
      m_kin_e(pts.size(), 0.0) {
  std::vector<double> weights_sqrt;
  for (const auto& w : m_weights) {
    weights_sqrt.push_back(sqrt(w));
  }

  auto weights_mat = negele::matrix::SquareMatrix::from_diag(weights_sqrt);
  negele::matrix::SquareMatrix mels_wo_weights(m_dim, mels);

  m_mels_with_weights = weights_mat * mels_wo_weights * weights_mat;

  std::vector<double> kin_diag;
  for (const auto& p : pts) {
    kin_diag.push_back(p * p * HBARC * HBARC / (2 * MASS));
  }
  m_kin_e = negele::matrix::SquareMatrix::from_diag(kin_diag);
}

void negele::matrix::MomentumSpaceMatrix::clear_mels() {
  m_mels_with_weights *= 0.0;
}
double negele::matrix::MomentumSpaceMatrix::max() const {
  return negele::matrix::abs(m_mels_with_weights).max();
}

void negele::matrix::MomentumSpaceMatrix::resize(
    const negele::matrix::MomentumSpaceMatrix& other) {
  m_dim = other.m_dim;
  m_pts = std::vector<double>(other.m_pts);
  m_weights = std::vector<double>(other.m_weights);

  m_mels_with_weights = negele::matrix::SquareMatrix(m_dim, 0.0);
  m_kin_e = negele::matrix::SquareMatrix(m_dim, 0.0);
}

bool negele::matrix::MomentumSpaceMatrix::same_size(
    const negele::matrix::MomentumSpaceMatrix& other) const {
  if (m_dim != other.m_dim) {
    return false;
  }
  for (int i = 0; i < m_dim; i++) {
    if (m_pts[i] != other.m_pts[i]) {
      return false;
    }
    if (m_weights[i] != other.m_weights[i]) {
      return false;
    }
  }
  if (m_mels_with_weights.size() != other.m_mels_with_weights.size()) {
    return false;
  }
  if (m_kin_e.size() != other.m_kin_e.size()) {
    return false;
  }
  return true;
}

negele::matrix::SquareMatrix
negele::matrix::MomentumSpaceMatrix::matrix_elements() const {
  std::vector<double> weights_inv_sqrt;
  for (const auto& w : m_weights) {
    weights_inv_sqrt.push_back(1.0 / sqrt(w));
  }
  auto weights_inv_mat =
      negele::matrix::SquareMatrix::from_diag(weights_inv_sqrt);
  return weights_inv_mat * m_mels_with_weights * weights_inv_mat;
}

std::vector<double> negele::matrix::mat_mul(const std::vector<double>& a,
                                            const std::vector<double>& b,
                                            int dim) {
  std::vector<double> new_matrix(dim * dim, 0.0);

  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      for (int k = 0; k < dim; k++) {
        new_matrix[i * dim + j] += a[i * dim + k] * b[k * dim + j];
      }
    }
  }

  return new_matrix;
}

negele::matrix::MomentumSpaceMatrix& negele::matrix::MomentumSpaceMatrix::
operator+=(const negele::matrix::MomentumSpaceMatrix& other) {
  m_mels_with_weights += other.m_mels_with_weights;
  return *this;
}

negele::matrix::MomentumSpaceMatrix& negele::matrix::MomentumSpaceMatrix::
operator+=(const double other) {
  for (int i = 0; i < m_dim * m_dim; i++) {
    m_mels_with_weights[i] += copysign(other, m_mels_with_weights[i]);
  }
  return *this;
}

negele::matrix::MomentumSpaceMatrix& negele::matrix::MomentumSpaceMatrix::
operator*=(const double factor) {
  m_mels_with_weights *= factor;
  return *this;
}

negele::matrix::MomentumSpaceMatrix negele::matrix::operator/(
    const negele::matrix::MomentumSpaceMatrix& a,
    const negele::matrix::MomentumSpaceMatrix& b) {
  negele::matrix::MomentumSpaceMatrix c(a);
  c.m_mels_with_weights /= b.m_mels_with_weights;
  return c;
}

negele::matrix::MomentumSpaceMatrix negele::matrix::abs(
    const negele::matrix::MomentumSpaceMatrix& a) {
  negele::matrix::MomentumSpaceMatrix c(a);
  c.m_mels_with_weights.abs();
  return c;
}

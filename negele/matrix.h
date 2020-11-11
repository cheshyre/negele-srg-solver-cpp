// Copyright 2020 Matthias Heinz
#ifndef NEGELE_MATRIX_H_
#define NEGELE_MATRIX_H_

#include <vector>

#include <boost/numeric/odeint.hpp>

namespace negele {
namespace matrix {

class SquareMatrix {
 private:
  int m_dim;
  std::vector<double> m_matrix;

 public:
  SquareMatrix(int dim, double value);

  SquareMatrix(int dim, const std::vector<double>& mels);

  static SquareMatrix from_diag(const std::vector<double>& diag_mels);

  SquareMatrix& operator+=(const SquareMatrix& other);
  SquareMatrix& operator-=(const SquareMatrix& other);
  SquareMatrix& operator*=(double factor);
  SquareMatrix& operator/=(double factor);
  SquareMatrix& operator/=(const SquareMatrix& other);
  SquareMatrix& abs();
  double max() const;
  const double& operator[](std::size_t index) const { return m_matrix[index]; }
  double& operator[](std::size_t index) { return m_matrix[index]; }
  std::size_t size() const { return m_dim * m_dim; }
  std::vector<double> eigenvalues() const;

  friend SquareMatrix operator*(const SquareMatrix& a, const SquareMatrix& b);
};

SquareMatrix operator+(const SquareMatrix& a, const SquareMatrix& b);
SquareMatrix operator-(const SquareMatrix& a, const SquareMatrix& b);
SquareMatrix operator*(double factor, const SquareMatrix& mat);
SquareMatrix operator*(const SquareMatrix& mat, double factor);
SquareMatrix operator/(const SquareMatrix& mat, double factor);
SquareMatrix operator/(const SquareMatrix& a, const SquareMatrix& b);
SquareMatrix abs(const SquareMatrix& a);

SquareMatrix operator*(const SquareMatrix& a, const SquareMatrix& b);
inline SquareMatrix commutator(const SquareMatrix& a, const SquareMatrix& b) {
  return a * b - b * a;
}

class MomentumSpaceMatrix {
 private:
  int m_dim;
  std::vector<double> m_pts;
  std::vector<double> m_weights;

  SquareMatrix m_mels_with_weights;
  SquareMatrix m_kin_e;

 public:
  friend MomentumSpaceMatrix operator/(const MomentumSpaceMatrix& a,
                                       const MomentumSpaceMatrix& b);
  friend MomentumSpaceMatrix abs(const MomentumSpaceMatrix& a);
  double max() const;
  /**
   * @brief Construct a new MomentumSpaceMatrix object.
   *
   * @param pts
   * @param weights
   * @param mels
   */
  MomentumSpaceMatrix(const std::vector<double>& pts,
                      const std::vector<double>& weights,
                      const std::vector<double>& mels);

  MomentumSpaceMatrix()
      : m_dim(0),
        m_pts(0, 0.0),
        m_weights(0, 0.0),
        m_mels_with_weights(0, 0.0),
        m_kin_e(0, 0.0) {}

  /**
   * @brief Set all the matrix elements to 0.
   *
   */
  void clear_mels();

  void resize(const MomentumSpaceMatrix& other);
  bool same_size(const MomentumSpaceMatrix& other) const;

  /**
   * @brief Get integration weights.
   *
   * @return std::vector<double>
   */
  std::vector<double> weights() const { return std::vector<double>(m_weights); }

  /**
   * @brief Get kinetic energy.
   *
   * @return SquareMatrix
   */
  SquareMatrix kinetic_energy() const { return m_kin_e; }

  /**
   * @brief Get the matrix elements with weights factored in.
   *
   * @return SquareMatrix
   */
  SquareMatrix matrix_elements_with_weights() const {
    return m_mels_with_weights;
  }

  /**
   * @brief Set the matrix elements (with weights included) in object.
   *
   * @param matrix
   */
  void set_matrix_elements(const SquareMatrix& matrix) {
    m_mels_with_weights = SquareMatrix(matrix);
  }

  /**
   * @brief Get the matrix elements (without integration weights).
   *
   * @return SquareMatrix
   */
  SquareMatrix matrix_elements() const;

  /**
   * @brief Add another matrix.
   *
   * @param other
   * @return MomentumSpaceMatrix&
   */
  MomentumSpaceMatrix& operator+=(const MomentumSpaceMatrix& other);

  /**
   * @brief Multiply by a factor.
   *
   * @param factor
   * @return MomentumSpaceMatrix&
   */
  MomentumSpaceMatrix& operator*=(const double factor);
  //   MomentumSpaceMatrix& operator/=(const MomentumSpaceMatrix& other);

  MomentumSpaceMatrix& operator+=(const double factor);
};

/**
 * @brief Evaluate matrix multiplication.
 *
 * @param a
 * @param b
 * @param dim
 * @return std::vector<double>
 */
std::vector<double> mat_mul(const std::vector<double>& a,
                            const std::vector<double>& b, int dim);

/**
 * @brief Add two matrices.
 *
 * @param a
 * @param b
 * @return MomentumSpaceMatrix
 */
inline MomentumSpaceMatrix operator+(const MomentumSpaceMatrix& a,
                                     const MomentumSpaceMatrix& b) {
  MomentumSpaceMatrix c(a);
  c += b;
  return c;
}

inline MomentumSpaceMatrix operator+(const MomentumSpaceMatrix& a,
                                     const double factor) {
  MomentumSpaceMatrix c(a);
  c += factor;
  return c;
}

inline MomentumSpaceMatrix operator+(const double factor,
                                     const MomentumSpaceMatrix& a) {
  MomentumSpaceMatrix c(a);
  c += factor;
  return c;
}

/**
 * @brief Multiply a matrix by a factor.
 *
 * @param a
 * @param b
 * @return MomentumSpaceMatrix
 */
inline MomentumSpaceMatrix operator*(const MomentumSpaceMatrix& a,
                                     const double b) {
  MomentumSpaceMatrix c(a);
  c *= b;
  return c;
}

/**
 * @brief Multiply a matrix by a factor.
 *
 * @param b
 * @param a
 * @return MomentumSpaceMatrix
 */
inline MomentumSpaceMatrix operator*(const double b,
                                     const MomentumSpaceMatrix& a) {
  MomentumSpaceMatrix c(a);
  c *= b;
  return c;
}

MomentumSpaceMatrix operator/(const MomentumSpaceMatrix& a,
                              const MomentumSpaceMatrix& b);
MomentumSpaceMatrix abs(const MomentumSpaceMatrix& a);

}  // namespace matrix
}  // namespace negele

namespace boost {
namespace numeric {
namespace odeint {
template <>
struct vector_space_norm_inf<negele::matrix::MomentumSpaceMatrix> {
  typedef double result_type;
  double operator()(const negele::matrix::MomentumSpaceMatrix& p) const {
    return p.max();
  }
};
}  // namespace odeint
}  // namespace numeric
}  // namespace boost

#endif  // NEGELE_MATRIX_H_

// Copyright 2020 Matthias Heinz
#ifndef NEGELE_POTENTIAL_H_
#define NEGELE_POTENTIAL_H_

namespace negele {
namespace potential {

class NegelePotential {
 private:
  double m_v1;
  double m_v2;
  double m_sig1;
  double m_sig2;

 public:
  /**
   * @brief Construct a new NegelePotential object.
   *
   * @param v1
   * @param v2
   * @param sig1
   * @param sig2
   */
  NegelePotential(double v1, double v2, double sig1, double sig2)
      : m_v1(v1), m_v2(v2), m_sig1(sig1), m_sig2(sig2) {}
  /**
   * @brief Evalute value of potential given momenta.
   *
   * @param p_in
   * @param p_out
   * @return double
   */
  double eval(double p_in, double p_out);
};

}  // namespace potential
}  // namespace negele

#endif  // NEGELE_POTENTIAL_H_

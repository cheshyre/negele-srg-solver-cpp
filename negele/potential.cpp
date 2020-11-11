// Copyright 2020 Matthias Heinz

#include "negele/potential.h"

#include <math.h>

double negele::potential::NegelePotential::eval(double p_in, double p_out) {
  return m_v1 / (2 * sqrt(2) * M_PI) *
             exp(-1 * (p_in - p_out) * (p_in - p_out) * m_sig1 * m_sig1 / 8.0) +
         m_v2 / (2 * sqrt(2) * M_PI) *
             exp(-1 * (p_in - p_out) * (p_in - p_out) * m_sig2 * m_sig2 / 8.0);
}

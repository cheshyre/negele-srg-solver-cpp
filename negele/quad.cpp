// Copyright 2020 Matthias Heinz
#include "negele/quad.h"

#include <math.h>
#include <tuple>
#include <vector>

std::tuple<std::vector<double>, std::vector<double>>
negele::quad::get_gauss_legendre_quadrature(int num_pts, double x_min,
                                            double x_max) {
  std::vector<double> x(num_pts, 0.0);
  std::vector<double> w(num_pts, 0.0);

  double eps = 10e-12;

  int m = (num_pts + 1) / 2;
  double xm = 0.5 * (x_max + x_min);
  double xl = 0.5 * (x_max - x_min);

  for (int i = 1; i <= m; i++) {
    double z = cos(M_PI * (i - 0.25) / (num_pts + 0.5));
    double z1 = 1.0;
    double pp = 0.0;

    while (fabs(z - z1) > eps) {
      double p1 = 1.0;
      double p2 = 0.0;

      for (int j = 1; j <= num_pts; j++) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }

      pp = num_pts * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }

    x[i - 1] = xm - xl * z;
    x[num_pts - i] = xm + xl * z;
    w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[num_pts - i] = w[i - 1];
  }

  return std::make_tuple(x, w);
}

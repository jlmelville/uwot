//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2018 James Melville
//
//  This file is part of UWOT
//
//  UWOT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  UWOT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with UWOT.  If not, see <http://www.gnu.org/licenses/>.
#include <cmath>
#include "gradient.h"

// UMAP

umap_gradient::umap_gradient(const double a, const double b, const double gamma) :
  a(a), b(b), a_b_m2(-2.0 * a * b), gamma_b_2(2.0 * gamma * b)
  { }

// The Python implementation does two pow calls: one for d^b, and one for
// d^(b-1) I have replaced the second call with (d^b) / d, which saves on the
// pow call at the cost of an extra multiply. There shouldn't be a risk of a
// divide by zero because I don't allow dist_squared to get lower than the
// machine limit double epsilon. Note that grad_rep has a (0.001 + dist_squared)
// term in its denominator, and I assume that the 0.001 is to prevent
// division-by-zero, so to be consistent with dealing with division by
// dist_squared I could add it here. But I haven't (to be more consistent with
// the original implementation).
// If this ever needs to change back:
// store a new const double field: const double b_m1 = b - 1.0
//  return (a_b_m2 * std::pow(dist_squared, b_m1)) / (a * std::pow(dist_squared, b) + 1.0);
const double umap_gradient::grad_attr(const double dist_squared) const {
  const double pd2b = std::pow(dist_squared, b);
  return (a_b_m2 * pd2b) / (dist_squared * (a * pd2b + 1.0));
}

const double umap_gradient::grad_rep(const double dist_squared) const {
  return gamma_b_2 / ((0.001 + dist_squared) * (a * std::pow(dist_squared, b) + 1.0));
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// an approximation to pow
double fastPrecisePow(double a, double b) {
  // calculate approximation with fraction of the exponent
  int e = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - e) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;

  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (e) {
    if (e & 1) {
      r *= a;
    }
    a *= a;
    e >>= 1;
  }

  return r * u.d;
}

// Approximate-Power UMAP: UMAP with an approximation to the slow pow calculation
apumap_gradient::apumap_gradient(const double a, const double b, const double gamma) :
  a(a), b(b), a_b_m2(-2.0 * a * b), gamma_b_2(2.0 * gamma * b)
  {  }

const double apumap_gradient::grad_attr(const double dist_squared) const {
  const double pd2b = fastPrecisePow(dist_squared, b);
  return (a_b_m2 * pd2b) / (dist_squared * (a * pd2b + 1.0));
}

const double apumap_gradient::grad_rep(const double dist_squared) const {
  return gamma_b_2 /
    ((0.001 + dist_squared) * (a * fastPrecisePow(dist_squared, b) + 1.0));
}

// t-UMAP

tumap_gradient::tumap_gradient() {}

const double tumap_gradient::grad_attr(const double dist_squared) const {
  return -2.0 / (dist_squared + 1.0);
}

const double tumap_gradient::grad_rep(const double dist_squared) const {
  return 2.0 / ((0.001 + dist_squared) * (dist_squared + 1.0));
}


// LargeVis

largevis_gradient::largevis_gradient(const double gamma) : 
  gamma_2(gamma * 2.0) {}

const double largevis_gradient::grad_attr(const double dist_squared) const {
  return -2.0 / (dist_squared + 1.0);
}

const double largevis_gradient::grad_rep(const double dist_squared) const {
  return gamma_2 / ((0.1 + dist_squared) * (dist_squared + 1.0));
}

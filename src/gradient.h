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

#ifndef UWOT_GRADIENT_H
#define UWOT_GRADIENT_H

// Class templated on the powfun function as suggested by Aaron Lun
template<double (*powfun)(double, double)>
class base_umap_gradient {
public:
  base_umap_gradient(const double a, const double b, const double gamma);
  const double grad_attr(const double dist_squared) const;
  const double grad_rep(const double dist_squared) const;
  static const constexpr double clamp_hi = 4.0;
  static const constexpr double clamp_lo = -4.0;
private:
  const double a;
  const double b;
  const double a_b_m2;
  const double gamma_b_2;
};

// UMAP
typedef base_umap_gradient<std::pow> umap_gradient;

// apUMAP: UMAP with an approximate power calculation
double fastPrecisePow(double, double);
typedef base_umap_gradient<fastPrecisePow> apumap_gradient;

// t-UMAP: the UMAP function with a = 1, and b = 1, which results in the Cauchy
// distribution as used in t-SNE. This massively simplifies the gradient,
// removing the pow calls, resulting in a noticeable speed increase (50% with
// MNIST), although the resulting embedding has a larger spread than the
// default. Also gamma is absent from this, because I believe it to be
// un-necessary in the UMAP cost function.
class tumap_gradient {
public:
  tumap_gradient();
  const double grad_attr(const double dist_squared) const;
  const double grad_rep(const double dist_squared) const;
  static const constexpr double clamp_hi = 4.0;
  static const constexpr double clamp_lo = -4.0;
};

class largevis_gradient {
public:
  largevis_gradient(const double gamma);
  const double grad_attr(const double dist_squared) const;
  const double grad_rep(const double dist_squared) const;
  static const constexpr double clamp_hi = 5.0;
  static const constexpr double clamp_lo = -5.0;
private:
  const double gamma_2;
};

#endif // UWOT_GRADIENT_H

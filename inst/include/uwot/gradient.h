// BSD 2-Clause License
//
// Copyright 2020 James Melville
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
//    this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// OF SUCH DAMAGE.

#ifndef UWOT_GRADIENT_H
#define UWOT_GRADIENT_H

#include <cmath>

namespace uwot {

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// an approximation to pow
inline auto fastPrecisePow(float a, float b) -> float {
  // calculate approximation with fraction of the exponent
  int e = static_cast<int>(b);
  union {
    double d;
    int x[2];
  } u = {a};
  u.x[1] = static_cast<int>((b - e) * (u.x[1] - 1072632447) + 1072632447);
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

  return static_cast<float>(r * u.d);
}

// Class templated on the powfun function as suggested by Aaron Lun
template <float (*powfun)(float, float)> class base_umap_gradient {
public:
  base_umap_gradient(float a, float b, float gamma)
      : a(a), b(b), a_b_m2(-2.0 * a * b), gamma_b_2(2.0 * gamma * b){};
  auto grad_attr(float dist_squared) const -> float {
    float pd2b = powfun(dist_squared, b);
    return (a_b_m2 * pd2b) / (dist_squared * (a * pd2b + 1.0));
  }
  auto grad_rep(float dist_squared) const -> float {
    return gamma_b_2 /
           ((0.001 + dist_squared) * (a * powfun(dist_squared, b) + 1.0));
  }

  static const constexpr float clamp_hi = 4.0;
  static const constexpr float clamp_lo = -4.0;

private:
  float a;
  float b;
  float a_b_m2;
  float gamma_b_2;
};

// UMAP using standard power function
using umap_gradient = base_umap_gradient<std::pow>;
// apUMAP: UMAP with an approximate power calculation
using apumap_gradient = base_umap_gradient<fastPrecisePow>;

// t-UMAP: the UMAP function with a = 1, and b = 1, which results in the Cauchy
// distribution as used in t-SNE. This massively simplifies the gradient,
// removing the pow calls, resulting in a noticeable speed increase (50% with
// MNIST), although the resulting embedding has a larger spread than the
// default. Also gamma is absent from this, because I believe it to be
// un-necessary in the UMAP cost function.
class tumap_gradient {
public:
  tumap_gradient() = default;
  auto grad_attr(float dist_squared) const -> float {
    return -2.0 / (dist_squared + 1.0);
  }
  auto grad_rep(float dist_squared) const -> float {
    return 2.0 / ((0.001 + dist_squared) * (dist_squared + 1.0));
  }
  static const constexpr float clamp_hi = 4.0;
  static const constexpr float clamp_lo = -4.0;
};

class largevis_gradient {
public:
  largevis_gradient(float gamma) : gamma_2(gamma * 2.0) {}
  auto grad_attr(float dist_squared) const -> float {
    return -2.0 / (dist_squared + 1.0);
  }
  auto grad_rep(float dist_squared) const -> float {
    return gamma_2 / ((0.1 + dist_squared) * (dist_squared + 1.0));
  }

  static const constexpr float clamp_hi = 5.0;
  static const constexpr float clamp_lo = -5.0;

private:
  float gamma_2;
};
} // namespace uwot

#endif // UWOT_GRADIENT_H

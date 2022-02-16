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
#include <limits>
#include <vector>

namespace uwot {

inline auto clamp(float v, float lo, float hi) -> float {
  float t = v < lo ? lo : v;
  return t > hi ? hi : t;
}

// return the squared euclidean distance between two points x[px] and y[py]
// also store the displacement between x[px] and y[py] in diffxy
// there is a small but noticeable performance improvement by doing so
// rather than recalculating it in the gradient step
inline auto d2diff(const std::vector<float> &x, std::size_t px,
                   const std::vector<float> &y, std::size_t py,
                   std::size_t ndim, float dist_eps, std::vector<float> &diffxy)
    -> float {
  float d2 = 0.0;
  for (std::size_t d = 0; d < ndim; d++) {
    float diff = x[px + d] - y[py + d];
    diffxy[d] = diff;
    d2 += diff * diff;
  }
  return (std::max)(dist_eps, d2);
}

// The gradient for the dth component of the displacement between two point,
// which for Euclidean distance in the output is invariably grad_coeff * (X - Y)
// Different methods clamp the magnitude of the gradient to different values
template <typename Gradient>
auto grad_d(const Gradient &gradient, const std::vector<float> &disp,
            std::size_t d, float grad_coeff) -> float {
  return gradient.clamp_grad(grad_coeff * disp[d]);
}

template <typename Gradient>
auto grad_attr(const Gradient &gradient,
               const std::vector<float> &head_embedding, std::size_t dj,
               const std::vector<float> &tail_embedding, std::size_t dk,
               std::size_t ndim, std::vector<float> &disp) -> float {
  static const float dist_eps = std::numeric_limits<float>::epsilon();
  float d2 =
      d2diff(head_embedding, dj, tail_embedding, dk, ndim, dist_eps, disp);
  return gradient.grad_attr(d2, dj, dk);
}

template <typename Gradient>
auto grad_rep(const Gradient &gradient,
              const std::vector<float> &head_embedding, std::size_t dj,
              const std::vector<float> &tail_embedding, std::size_t dk,
              std::size_t ndim, std::vector<float> &disp) -> float {
  static const float dist_eps = std::numeric_limits<float>::epsilon();
  float d2 =
      d2diff(head_embedding, dj, tail_embedding, dk, ndim, dist_eps, disp);
  return gradient.grad_rep(d2, dj, dk);
}

// https://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
// an approximation to pow
inline auto fastPrecisePow(float a, float b) -> float {
  // calculate approximation with fraction of the exponent
  int e = static_cast<int>(b);
  union {
    double d;
    int x[2];
  } u = {a};
  u.x[1] = static_cast<int>((b - static_cast<double>(e)) *
                                static_cast<double>(u.x[1] - 1072632447) +
                            1072632447.0);
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
  // Compared to the UMAP Python implementation, instead of doing d2^(b-1)
  // we can save a power calculation by using d2^b / d2
  auto grad_attr(float d2, std::size_t, std::size_t) const -> float {
    float pd2b = powfun(d2, b);
    return (a_b_m2 * pd2b) / (d2 * (a * pd2b + 1.0));
  }
  auto grad_rep(float d2, std::size_t, std::size_t) const -> float {
    return gamma_b_2 / ((0.001 + d2) * (a * powfun(d2, b) + 1.0));
  }
  auto clamp_grad(float grad_d) const -> float {
    return clamp(grad_d, clamp_lo, clamp_hi);
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
  auto grad_attr(float d2, std::size_t, std::size_t) const -> float {
    return -2.0 / (d2 + 1.0);
  }
  auto grad_rep(float d2, std::size_t, std::size_t) const -> float {
    return 2.0 / ((0.001 + d2) * (d2 + 1.0));
  }
  auto clamp_grad(float grad_d) const -> float {
    return clamp(grad_d, clamp_lo, clamp_hi);
  }
  static const constexpr float clamp_hi = 4.0;
  static const constexpr float clamp_lo = -4.0;
};

// UMAP where a varies for each observation
class umapai_gradient {
public:
  umapai_gradient(const std::vector<float> &ai, float b, std::size_t ndim)
      : ai(ai), b(b), ndim(ndim), b_m2(-2.0 * b), b_2(2.0 * b) {}

  auto clamp_grad(float grad_d) const -> float {
    return clamp(grad_d, clamp_lo, clamp_hi);
  }
  static const constexpr float clamp_hi = 4.0;
  static const constexpr float clamp_lo = -4.0;
  auto grad_attr(float d2, std::size_t i, std::size_t j) const -> float {
    auto a = ai[i / ndim] * ai[j / ndim];
    float pd2b = std::pow(d2, b);
    return (a * b_m2 * pd2b) / (d2 * (a * pd2b + 1.0));
  }

  auto grad_rep(float d2, std::size_t i, std::size_t j) const -> float {
    auto a = ai[i / ndim] * ai[j / ndim];
    return b_2 / ((0.001 + d2) * (a * std::pow(d2, b) + 1.0));
  }

private:
  std::vector<float> ai;
  float b;
  std::size_t ndim;
  float b_m2;
  float b_2;
};

class umapai2_gradient {
public:
  umapai2_gradient(const std::vector<float> &ai, const std::vector<float> &aj,
                   float b, std::size_t ndim)
      : ai(ai), aj(aj), b(b), ndim(ndim), b_m2(-2.0 * b), b_2(2.0 * b) {}

  auto clamp_grad(float grad_d) const -> float {
    return clamp(grad_d, clamp_lo, clamp_hi);
  }
  static const constexpr float clamp_hi = 4.0;
  static const constexpr float clamp_lo = -4.0;
  auto grad_attr(float d2, std::size_t i, std::size_t j) const -> float {
    auto a = ai[i / ndim] * aj[j / ndim];
    float pd2b = std::pow(d2, b);
    return (a * b_m2 * pd2b) / (d2 * (a * pd2b + 1.0));
  }

  auto grad_rep(float d2, std::size_t i, std::size_t j) const -> float {
    auto a = ai[i / ndim] * aj[j / ndim];
    return b_2 / ((0.001 + d2) * (a * std::pow(d2, b) + 1.0));
  }

private:
  std::vector<float> ai;
  std::vector<float> aj;
  float b;
  std::size_t ndim;
  float b_m2;
  float b_2;
};

class largevis_gradient {
public:
  largevis_gradient(float gamma) : gamma_2(gamma * 2.0) {}
  auto grad_attr(float d2, std::size_t, std::size_t) const -> float {
    return -2.0 / (d2 + 1.0);
  }
  auto grad_rep(float d2, std::size_t, std::size_t) const -> float {
    return gamma_2 / ((0.1 + d2) * (d2 + 1.0));
  }
  auto clamp_grad(float grad_d) const -> float {
    return clamp(grad_d, clamp_lo, clamp_hi);
  }
  static const constexpr float clamp_hi = 5.0;
  static const constexpr float clamp_lo = -5.0;

private:
  float gamma_2;
};

class pacmap_gradient {
public:
  float a;
  float b;
  float ab2m;
  float b1;

  // near: a 1-3, b = 10; mid: a = 3-1000, b = 10,000
  pacmap_gradient(float a, float b)
      : a(1.0), b(10.0), ab2m(-2.0 * a * b), b1(b + 1.0) {}
  auto grad_attr(float d2, std::size_t, std::size_t) const -> float {
    return ab2m / ((b1 + d2) * (b1 + d2));
  }
  auto grad_rep(float d2, std::size_t, std::size_t) const -> float {
    return 2.0 / ((2.0 + d2) * (2.0 + d2));
  }
  inline auto clamp_grad(float grad_d) const -> float { return grad_d; }
};
} // namespace uwot

#endif // UWOT_GRADIENT_H

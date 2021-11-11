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

#ifndef UWOT_OPTIMIZE_H
#define UWOT_OPTIMIZE_H

#include <vector>

namespace uwot {

float linear_decay(double val, std::size_t epoch, std::size_t n_epochs) {
  return val *
         (1.0 - (static_cast<float>(epoch) / static_cast<float>(n_epochs)));
}

float linear_grow(double val, std::size_t epoch, std::size_t n_epochs) {
  return val * (static_cast<float>(epoch) / static_cast<float>(n_epochs));
}

struct Sgd {
  float initial_alpha;
  float alpha;

  Sgd(float alpha) : initial_alpha(alpha), alpha(alpha){};

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    v[i] += alpha * grad[i];
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha = linear_decay(initial_alpha, epoch, n_epochs);
  }
};

struct Adam {
  float initial_alpha;
  float alpha;
  float beta1;
  float beta2;

  float beta11; // 1 - beta1
  float beta1t; // beta1 ^ t

  float beta21; // 1 - beta2
  float beta2t; // beta2 ^ t

  float eps;

  // rather than calculate the debiased values for m and v (mhat and vhat)
  // directly, fold them into a scaling factor for alpha as described at the
  // end of the first paragraph of section 2 of the Adam paper
  // technically you also need to rescale eps (given as eps-hat in the paper)
  float epsc;     // scaled eps
  float ad_scale; // scaled alpha

  std::vector<float> mt;
  std::vector<float> vt;

  Adam(float alpha, float beta1, float beta2, float eps, std::size_t vec_size)
      : initial_alpha(alpha), alpha(alpha), beta1(beta1), beta2(beta2),
        beta11(1.0 - beta1), beta1t(beta1), beta21(1.0 - beta2), beta2t(beta2),
        eps(eps), epsc(eps * sqrt(beta21)),
        ad_scale(alpha * sqrt(beta21) / beta11), mt(vec_size), vt(vec_size) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    // this takes advantage of updating in-place to give a more compact version
    // of mt[i] = beta1 * mt[i] + beta11 * grad[i] etc.
    vt[i] += beta21 * (grad[i] * grad[i] - vt[i]);
    mt[i] += beta11 * (grad[i] - mt[i]);

    // ad_scale and epsc handle the debiasing
    v[i] += ad_scale * mt[i] / (sqrt(vt[i]) + epsc);
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha = linear_decay(initial_alpha, epoch, n_epochs);

    // update debiasing factors
    beta1t *= beta1;
    beta2t *= beta2;
    float sqrt_b2t1 = sqrt(1.0 - beta2t);

    // rescale alpha and eps to take account of debiasing
    ad_scale = alpha * sqrt_b2t1 / (1.0 - beta1t);
    epsc = sqrt_b2t1 * eps;
  }
};

} // namespace uwot

#endif // UWOT_OPTIMIZE_H

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

#include <iostream>
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

struct ParamUpdate {
  virtual float update(float val, std::size_t, std::size_t) = 0;
};

struct ParamFixed : ParamUpdate {
  virtual float update(float val, std::size_t, std::size_t) override {
    return val;
  }
};

struct ParamLinearGrow : ParamUpdate {
  virtual float update(float val, std::size_t epoch,
                       std::size_t n_epochs) override {
    return val * (static_cast<float>(epoch) / static_cast<float>(n_epochs));
  }
};

struct ParamLinearDecay : ParamUpdate {
  virtual float update(float val, std::size_t epoch,
                       std::size_t n_epochs) override {
    return val *
           (1.0 - (static_cast<float>(epoch) / static_cast<float>(n_epochs)));
  }
};

struct ParamDemon : ParamUpdate {
  virtual float update(float val, std::size_t epoch, std::size_t n_epochs) {
    const float tau =
        (1.0 - (static_cast<float>(epoch)) / static_cast<float>(n_epochs));
    return val * (1.0 - tau) / (1.0 - val * tau);
  }
};

struct ParamSlowStart : ParamUpdate {
  std::size_t n_start_epochs;
  std::unique_ptr<ParamUpdate> param_update;

  ParamSlowStart(std::size_t n_start_epochs, ParamUpdate *param_update)
      : n_start_epochs(n_start_epochs), param_update(param_update) {}

  virtual float update(float val, std::size_t epoch,
                       std::size_t n_epochs) override {
    float updated = param_update->update(val, epoch, n_epochs);
    return std::min(updated, (static_cast<float>(epoch) /
                              static_cast<float>(n_start_epochs)) *
                                 updated);
  }
};

struct MomentumSgd {
  float initial_alpha;
  float alpha;
  float initial_mu;
  float mu;
  float mu1;

  std::unique_ptr<uwot::ParamUpdate> alpha_param;
  std::unique_ptr<uwot::ParamUpdate> mu_param;

  std::vector<float> up_old;

  MomentumSgd(float alpha, float mu, std::size_t vec_len,
              uwot::ParamUpdate *alpha_update, uwot::ParamUpdate *mu_update)
      : initial_alpha(alpha), alpha(alpha), initial_mu(mu), mu(mu),
        mu1(1.0 - mu), alpha_param(std::move(alpha_update)),
        mu_param(std::move(mu_update)), up_old(vec_len) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    float up = mu1 * grad[i] + mu * up_old[i];
    v[i] += alpha * up;
    up_old[i] = up;
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha = alpha_param->update(initial_alpha, epoch, n_epochs);
    mu = mu_param->update(initial_mu, epoch, n_epochs);
    mu1 = 1.0 - mu;
  }
};

struct Adam {
  float initial_alpha;
  float alpha;

  float initial_beta1;
  float beta1;
  float beta11;
  float beta1t;
  float beta1t1;

  float beta2;
  float beta21;
  float beta2t;
  float beta2t1;

  float eps;

  std::vector<float> mt;
  std::vector<float> vt;

  std::unique_ptr<uwot::ParamUpdate> alpha_param;
  std::unique_ptr<uwot::ParamUpdate> beta1_param;

  Adam(float alpha, float beta1, float beta2, float eps, std::size_t vec_size,
       uwot::ParamUpdate *alpha_update, uwot::ParamUpdate *beta_update)
      : initial_alpha(alpha), alpha(alpha), initial_beta1(beta1), beta1(beta1),
        beta11(1.0 - beta1), beta1t(beta1), beta1t1(1.0 - beta1t), beta2(beta2),
        beta21(1.0 - beta2), beta2t(beta2), beta2t1(1.0 - beta2t), eps(eps),
        mt(vec_size), vt(vec_size), alpha_param(std::move(alpha_update)),
        beta1_param(std::move(beta_update)) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    float mb = beta1 * mt[i] + beta11 * grad[i];
    float vb = beta2 * vt[i] + beta21 * grad[i] * grad[i];

    float mc = mb / beta1t1;
    float vc = vb / beta2t1;

    float up = mc / (sqrt(vc) + eps);
    v[i] += alpha * up;
    mt[i] = mb;
    vt[i] = vb;
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha = alpha_param->update(initial_alpha, epoch, n_epochs);
    beta1 = beta1_param->update(initial_beta1, epoch, n_epochs);

    beta1t *= beta1;
    beta1t1 = 1.0 - beta1t;
    beta2t *= beta2;
    beta2t1 = 1.0 - beta2t;
  }
};

} // namespace uwot

#endif // UWOT_OPTIMIZE_H

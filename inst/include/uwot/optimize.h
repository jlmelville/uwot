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
    const float btau1 = val * (1.0 - (static_cast<float>(epoch)) /
                                         static_cast<float>(n_epochs));
    return btau1 / (1.0 - val + btau1);
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

struct Param {
  float initial_value;
  float value;
  std::unique_ptr<ParamUpdate> param_update;

  Param(float initial_value, ParamUpdate *param_update)
      : initial_value(initial_value), value(initial_value),
        param_update(std::move(param_update)) {}
  Param(Param &&other)
      : initial_value(other.initial_value), value(other.value),
        param_update(std::move(other.param_update)) {}

  void update(std::size_t epoch, std::size_t n_epochs) {
    value = param_update->update(initial_value, epoch, n_epochs);
  }
};

struct MomentumSgd {
  uwot::Param alpha_param;
  uwot::Param mu_param;

  float mu1;

  std::vector<float> up_old;

  MomentumSgd(uwot::Param &alpha_param, uwot::Param &mu_param,
              std::size_t vec_len)
      : alpha_param(std::move(alpha_param)), mu_param(std::move(mu_param)),
        mu1(1.0 - this->mu_param.value), up_old(vec_len) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    float up = mu1 * grad[i] + mu_param.value * up_old[i];
    v[i] += alpha_param.value * up;
    up_old[i] = up;
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha_param.update(epoch, n_epochs);
    mu_param.update(epoch, n_epochs);
    mu1 = 1.0 - mu_param.value;
    // std::cout << epoch << ": alpha = " << alpha_param.value
    //           << " beta = " << mu_param.value << std::endl;
  }
};

struct Adam {
  uwot::Param alpha_param;
  uwot::Param beta1_param;
  uwot::Param beta2_param;

  float beta11;
  float beta1t;
  float beta1t1;

  float beta21;
  float beta2t;
  float beta2t1;

  float eps;

  std::vector<float> mt;
  std::vector<float> vt;

  Adam(uwot::Param &alpha_param, uwot::Param &beta1_param,
       uwot::Param &beta2_param, float eps, std::size_t vec_size)
      : alpha_param(std::move(alpha_param)),
        beta1_param(std::move(beta1_param)),
        beta2_param(std::move(beta2_param)),
        beta11(1.0 - this->beta1_param.value), beta1t(this->beta1_param.value),
        beta1t1(1.0 - beta1t), beta21(1.0 - this->beta2_param.value),
        beta2t(this->beta2_param.value), beta2t1(1.0 - beta2t), eps(eps),
        mt(vec_size), vt(vec_size) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    float mb = beta1_param.value * mt[i] + beta11 * grad[i];
    float vb = beta2_param.value * vt[i] + beta21 * grad[i] * grad[i];

    float mc = mb / beta1t1;
    float vc = vb / beta2t1;

    float up = mc / (sqrt(vc) + eps);
    v[i] += alpha_param.value * up;
    mt[i] = mb;
    vt[i] = vb;
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha_param.update(epoch, n_epochs);
    beta1_param.update(epoch, n_epochs);

    beta1t *= beta1_param.value;
    beta1t1 = 1.0 - beta1t;
    beta2t *= beta2_param.value;
    beta2t1 = 1.0 - beta2t;

    // std::cout << epoch << ": alpha = " << alpha_param.value
    //           << " beta1 = " << beta1_param.value
    //           << " beta2 = " << beta2_param.value << std::endl;
  }
};

// Quasi-Hyperbolic Momentum (Ma & Yarats 2019)
// nu = 0: SGD
// nu = beta: NAG
// nu = 1: CM
// m_t = beta * m_t-1 + (1 - beta) * s_t-1 (CM)
// q_t = nu * m_t + (1 - nu) * s_t-1
// q_t = nu * beta * m_t-1 + (1 - nu * beta) * s_t-1 (expanded version)
struct Qhm {
  uwot::Param alpha_param;
  uwot::Param beta_param;
  uwot::Param nu_param;

  float beta1;
  float nu1;

  float bt; // beta products
  float nbt1; // (1 - nu * bt...b3b2b1)

  std::vector<float> mt;

  Qhm(uwot::Param &alpha_param, uwot::Param &beta_param, uwot::Param &nu_param,
      std::size_t vec_size)
      : alpha_param(std::move(alpha_param)), beta_param(std::move(beta_param)),
        nu_param(std::move(nu_param)), beta1(1.0 - this->beta_param.value),
        nu1(1.0 - this->nu_param.value),
        bt(this->beta_param.value),
        nbt1(1.0 - (this->nu_param.value * this->bt)), mt(vec_size) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    // float mold = mt[0];
    mt[i] = beta_param.value * mt[i] + beta1 * grad[i];
    float qb = nu_param.value * mt[i] + nu1 * grad[i];

    float up = qb / nbt1;
    v[i] += alpha_param.value * up;

    // if (i == 0) {
    //   std::cout << mold << " " << mt[0] << " " << qb << " " << up << " " << v[0] << std::endl;
    // }
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    // float nu_old = nu_param.value;
    
    alpha_param.update(epoch, n_epochs);
    beta_param.update(epoch, n_epochs);
    nu_param.update(epoch, n_epochs);

    beta1 = 1.0 - beta_param.value;
    nu1 = 1.0 - nu_param.value;

    bt *= beta_param.value;
    nbt1 = 1.0 - nu_param.value * bt;

    // std::cout << epoch << ": alpha = " << alpha_param.value
    //           << " beta = " << beta_param.value
    //           << " nu = " << nu_param.value << std::endl;
  }
};

// Quasi-Quasi-Hyperbolic Momentum
// like QHM but nu is scaled differently so its meaning is independent of beta
// nu = 1: CM
// nu = 2: NAG
// nu -> Inf: SGD (but just set nu=1, beta=0)
// nu < 1: places lower weight on most recent gradient
// nu is the number of averaging steps: e.g. once for CM, twice for NAG etc.
// qq_t = beta ^ nu m_t-1 + (1 - beta ^ nu) * s_t-1
// (easier to write in terms of the expanded CM step than the two-step QHM form)
struct Qqhm {
  uwot::Param alpha_param;
  uwot::Param beta_param;
  uwot::Param nu_param;

  // momentum weighting
  float beta1; // 1 - beta

  // iterative averaging
  float bn;  // beta ^ nu
  float bn1; // 1 - beta ^ nu

  // debiasing denominator
  float bt;  // bt...b3b2b1
  float bnt1; // 1.0 - bt^nu...b3b2b1
  std::vector<float> mt;

  Qqhm(uwot::Param &alpha_param, uwot::Param &beta_param, uwot::Param &nu_param,
       std::size_t vec_size)
      : alpha_param(std::move(alpha_param)), beta_param(std::move(beta_param)),
        nu_param(std::move(nu_param)), beta1(1.0 - this->beta_param.value),
        bn(std::pow(this->beta_param.value, this->nu_param.value)),
        bn1(1.0 - this->bn), 
        bt(this->beta_param.value), 
        bnt1(1.0 - this->bn),
        mt(vec_size) {}

  void update(std::vector<float> &v, std::vector<float> &grad, std::size_t i) {
    // float mold = mt[0];
    float qqb = bn * mt[i] + bn1 * grad[i];
  
    // debias
    float up = qqb / bnt1;
    
    // do the update
    v[i] += alpha_param.value * up;

    // momentum update for next step
    mt[i] = beta_param.value * mt[i] + beta1 * grad[i];
    
    // if (i == 0) {
    //   std::cout << mold << " " << mt[0] << " " << qqb << " " << up << " " << v[0] << std::endl;
    // }
  }

  void epoch_end(std::size_t epoch, std::size_t n_epochs) {
    alpha_param.update(epoch, n_epochs);
    beta_param.update(epoch, n_epochs);
    nu_param.update(epoch, n_epochs);

    beta1 = 1.0 - beta_param.value;
    bn = std::pow(beta_param.value, nu_param.value);
    bn1 = 1.0 - bn;

    bt *= beta_param.value; // ...b3b2b1
    // 1 - (bt...b3b2b1 x bt^(nu - 1)) = 1 - (bt^nu...b3b2b1)
    bnt1 = 1.0 - bt * std::pow(beta_param.value, nu_param.value - 1);
    

    // std::cout << epoch << ": alpha = " << alpha_param.value
    //           << " beta = " << beta_param.value
    //           << " nu = " << nu_param.value << std::endl;
  }
};

} // namespace uwot

#endif // UWOT_OPTIMIZE_H

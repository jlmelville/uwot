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

#include <vector>

#include "uwot/coords.h"
#include "uwot/epoch.h"
#include "uwot/gradient.h"
#include "uwot/optimize.h"
#include "uwot/sampler.h"

#include "rng.h"
#include "rparallel.h"
#include "rprogress.h"

using namespace Rcpp;

template <typename T>
auto lget(List list, const std::string &name, T default_value) -> T {
  auto key = name.c_str();
  if (!list.containsElementNamed(key)) {
    return default_value;
  } else {
    return list[key];
  }
}

// Template class specialization to handle different rng/batch combinations
template <bool DoBatch = true> struct BatchRngFactory {
  using PcgFactoryType = batch_pcg_factory;
  using TauFactoryType = batch_tau_factory;
};
template <> struct BatchRngFactory<false> {
  using PcgFactoryType = pcg_factory;
  using TauFactoryType = tau_factory;
};

struct UmapFactory {
  bool move_other;
  bool pcg_rand;
  std::vector<float> &head_embedding;
  std::vector<float> &tail_embedding;
  const std::vector<unsigned int> &positive_head;
  const std::vector<unsigned int> &positive_tail;
  const std::vector<unsigned int> &positive_ptr;
  unsigned int n_epochs;
  unsigned int n_head_vertices;
  unsigned int n_tail_vertices;
  const std::vector<float> &epochs_per_sample;
  float initial_alpha;
  List opt_args;
  float negative_sample_rate;
  bool batch;
  std::size_t n_threads;
  std::size_t grain_size;
  bool verbose;

  UmapFactory(bool move_other, bool pcg_rand,
              std::vector<float> &head_embedding,
              std::vector<float> &tail_embedding,
              const std::vector<unsigned int> &positive_head,
              const std::vector<unsigned int> &positive_tail,
              const std::vector<unsigned int> &positive_ptr,
              unsigned int n_epochs, unsigned int n_head_vertices,
              unsigned int n_tail_vertices,
              const std::vector<float> &epochs_per_sample, float initial_alpha,
              List opt_args, float negative_sample_rate, bool batch,
              std::size_t n_threads, std::size_t grain_size, bool verbose)
      : move_other(move_other), pcg_rand(pcg_rand),
        head_embedding(head_embedding), tail_embedding(tail_embedding),
        positive_head(positive_head), positive_tail(positive_tail),
        positive_ptr(positive_ptr), n_epochs(n_epochs),
        n_head_vertices(n_head_vertices), n_tail_vertices(n_tail_vertices),
        epochs_per_sample(epochs_per_sample), initial_alpha(initial_alpha),
        opt_args(opt_args), negative_sample_rate(negative_sample_rate),
        batch(batch), n_threads(n_threads), grain_size(grain_size),
        verbose(verbose) {}

  template <typename Gradient> void create(const Gradient &gradient) {
    uwot::Sampler sampler(epochs_per_sample, negative_sample_rate);

    if (move_other) {
      create_impl<true>(gradient, sampler, pcg_rand, batch);
    } else {
      create_impl<false>(gradient, sampler, pcg_rand, batch);
    }
  }

  template <bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, uwot::Sampler &sampler,
                   bool pcg_rand, bool batch) {
    if (batch) {
      create_impl<BatchRngFactory<true>, DoMove>(gradient, sampler, pcg_rand,
                                                 batch);
    } else {
      create_impl<BatchRngFactory<false>, DoMove>(gradient, sampler, pcg_rand,
                                                  batch);
    }
  }

  template <typename BatchRngFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, uwot::Sampler &sampler,
                   bool pcg_rand, bool batch) {
    if (pcg_rand) {
      create_impl<typename BatchRngFactory::PcgFactoryType, DoMove>(
          gradient, sampler, batch);
    } else {
      create_impl<typename BatchRngFactory::TauFactoryType, DoMove>(
          gradient, sampler, batch);
    }
  }

  auto create_param_update(const std::string &param_update_name)
      -> uwot::ParamUpdate * {
    if (param_update_name == "constant") {
      return new uwot::ParamFixed();
    } else if (param_update_name == "linear_grow") {
      return new uwot::ParamLinearGrow();
    } else if (param_update_name == "demon") {
      return new uwot::ParamDemon();
    } else if (param_update_name == "linear_decay") {
      return new uwot::ParamLinearDecay();
    } else {
      stop("Unknown param update");
    }
  }

  auto create_param_slow_start(uwot::ParamUpdate *param_update,
                               float slow_start_frac) -> uwot::ParamUpdate * {
    std::size_t start_epochs =
        std::ceil(slow_start_frac * static_cast<float>(n_epochs));
    Rcerr << " slow start for " << start_epochs << " epochs";
    return new uwot::ParamSlowStart(start_epochs, param_update);
  }

  auto create_param_update(List opt_args, const std::string &prefix,
                           const std::string &default_update)
      -> uwot::ParamUpdate * {
    std::string update = lget(opt_args, prefix + "_update", default_update);
    Rcerr << " " << prefix << " " << update;

    uwot::ParamUpdate *param_update = create_param_update(update);

    float slow_start = lget(opt_args, prefix + "_slow_start", 0.0);
    if (slow_start < 0.0 || slow_start > 1.0) {
      stop("slow start must be between 0-1");
    }
    if (slow_start > 0.0) {
      param_update = create_param_slow_start(param_update, slow_start);
    }

    return param_update;
  }

  auto create_adam(List opt_args) -> uwot::Adam {
    float alpha = lget(opt_args, "alpha", 1.0);
    float beta1 = lget(opt_args, "beta1", 0.9);
    float beta2 = lget(opt_args, "beta2", 0.999);
    float eps = lget(opt_args, "eps", 1e-8);

    Rcerr << "Optimizing with Adam alpha = " << alpha << " beta1 = " << beta1
          << " beta2 = " << beta2 << " eps = " << eps;

    uwot::ParamUpdate *alpha_update =
        create_param_update(opt_args, "alpha", "linear_decay");
    uwot::ParamUpdate *beta_update =
        create_param_update(opt_args, "beta1", "constant");
    Rcerr << std::endl;

    return uwot::Adam(alpha, beta1, beta2, eps, head_embedding.size(),
                      alpha_update, beta_update);
  }

  auto create_msgd(List opt_args) -> uwot::MomentumSgd {
    float alpha = lget(opt_args, "alpha", 1.0);
    float beta = lget(opt_args, "beta", 0.0);
    Rcerr << "Optimizing with momentum sgd alpha = " << alpha
          << " beta = " << beta;
    uwot::ParamUpdate *alpha_update =
        create_param_update(opt_args, "alpha", "linear_decay");
    uwot::ParamUpdate *beta_update =
        create_param_update(opt_args, "beta", "constant");
    Rcerr << std::endl;

    return uwot::MomentumSgd(alpha, beta, head_embedding.size(), alpha_update,
                             beta_update);
  }

  template <typename RandFactory, bool DoMove, typename Gradient>
  void create_impl(const Gradient &gradient, uwot::Sampler &sampler,
                   bool batch) {
    const std::size_t ndim = head_embedding.size() / n_head_vertices;

    if (batch) {
      std::string opt_name = opt_args[0];
      if (opt_name == "adam") {
        auto opt = create_adam(opt_args);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, sampler, opt, batch);
      } else if (opt_name == "msgd") {
        auto opt = create_msgd(opt_args);
        create_impl_batch_opt<decltype(opt), RandFactory, DoMove, Gradient>(
            gradient, sampler, opt, batch);
      } else {
        stop("Unknown optimization method");
      }
    } else {
      uwot::InPlaceUpdate<DoMove> update(head_embedding, tail_embedding,
                                         initial_alpha);
      uwot::EdgeWorker<Gradient, decltype(update), RandFactory> worker(
          gradient, update, positive_head, positive_tail, sampler, ndim,
          n_tail_vertices, n_threads);
      create_impl(worker, gradient, sampler);
    }
  }

  template <typename Opt, typename RandFactory, bool DoMove, typename Gradient>
  void create_impl_batch_opt(const Gradient &gradient, uwot::Sampler &sampler,
                             Opt &opt, bool batch) {
    uwot::BatchUpdate<DoMove, decltype(opt)> update(head_embedding,
                                                    tail_embedding, opt);
    const std::size_t ndim = head_embedding.size() / n_head_vertices;
    uwot::NodeWorker<Gradient, decltype(update), RandFactory> worker(
        gradient, update, positive_head, positive_tail, positive_ptr, sampler,
        ndim, n_tail_vertices);
    create_impl(worker, gradient, sampler);
  }

  template <typename Worker, typename Gradient>
  void create_impl(Worker &worker, const Gradient &gradient,
                   uwot::Sampler &sampler) {

    RProgress progress(n_epochs, verbose);
    if (n_threads > 0) {
      RParallel parallel(n_threads, grain_size);
      create_impl(worker, gradient, sampler, progress, parallel);
    } else {
      RSerial serial;
      create_impl(worker, gradient, sampler, progress, serial);
    }
  }

  template <typename Worker, typename Gradient, typename Progress,
            typename Parallel>
  void create_impl(Worker &worker, const Gradient &gradient,
                   uwot::Sampler &sampler, Progress &progress,
                   Parallel &parallel) {
    uwot::optimize_layout(worker, progress, n_epochs, parallel);
  }
};

auto r_to_coords(NumericMatrix head_embedding,
                 Nullable<NumericMatrix> tail_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  if (tail_embedding.isNull()) {
    return uwot::Coords(head_vec);
  } else {
    auto tail_vec = as<std::vector<float>>(tail_embedding);
    return uwot::Coords(head_vec, tail_vec);
  }
}

auto r_to_coords(NumericMatrix head_embedding) -> uwot::Coords {
  auto head_vec = as<std::vector<float>>(head_embedding);
  return uwot::Coords(head_vec);
}

void validate_args(List method_args,
                   const std::vector<std::string> &arg_names) {
  for (auto &arg_name : arg_names) {
    if (!method_args.containsElementNamed(arg_name.c_str())) {
      stop("Missing embedding method argument: " + arg_name);
    }
  }
}

void create_umap(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"a", "b", "gamma", "approx_pow"};
  validate_args(method_args, arg_names);

  float a = method_args["a"];
  float b = method_args["b"];
  float gamma = method_args["gamma"];
  bool approx_pow = method_args["approx_pow"];
  if (approx_pow) {
    const uwot::apumap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  } else {
    const uwot::umap_gradient gradient(a, b, gamma);
    umap_factory.create(gradient);
  }
}

void create_tumap(UmapFactory &umap_factory, List) {
  const uwot::tumap_gradient gradient;
  umap_factory.create(gradient);
}

void create_largevis(UmapFactory &umap_factory, List method_args) {
  std::vector<std::string> arg_names = {"gamma"};
  validate_args(method_args, arg_names);

  float gamma = method_args["gamma"];
  const uwot::largevis_gradient gradient(gamma);
  umap_factory.create(gradient);
}

// [[Rcpp::export]]
NumericMatrix optimize_layout_r(
    NumericMatrix head_embedding, Nullable<NumericMatrix> tail_embedding,
    const std::vector<unsigned int> positive_head,
    const std::vector<unsigned int> positive_tail,
    const std::vector<unsigned int> positive_ptr, unsigned int n_epochs,
    unsigned int n_head_vertices, unsigned int n_tail_vertices,
    const std::vector<float> epochs_per_sample, const std::string &method,
    List method_args, float initial_alpha, List opt_args,
    float negative_sample_rate, bool pcg_rand = true, bool batch = false,
    std::size_t n_threads = 0, std::size_t grain_size = 1,
    bool move_other = true, bool verbose = false) {

  auto coords = r_to_coords(head_embedding, tail_embedding);

  UmapFactory umap_factory(move_other, pcg_rand, coords.get_head_embedding(),
                           coords.get_tail_embedding(), positive_head,
                           positive_tail, positive_ptr, n_epochs,
                           n_head_vertices, n_tail_vertices, epochs_per_sample,
                           initial_alpha, opt_args, negative_sample_rate, batch,
                           n_threads, grain_size, verbose);

  if (method == "umap") {
    create_umap(umap_factory, method_args);
  } else if (method == "tumap") {
    create_tumap(umap_factory, method_args);
  } else if (method == "largevis") {
    create_largevis(umap_factory, method_args);
  } else {
    stop("Unknown method: '" + method + "'");
  }

  return NumericMatrix(head_embedding.nrow(), head_embedding.ncol(),
                       coords.get_head_embedding().begin());
}

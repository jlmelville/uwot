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

#ifndef UWOT_SAMPLER_H
#define UWOT_SAMPLER_H

#include <vector>

namespace uwot {

// Weighted edge sampler
class Sampler {
public:
  std::size_t epoch;

  Sampler(const std::vector<float> &epochs_per_sample,
          float negative_sample_rate)
      : epoch(0), epochs_per_sample(epochs_per_sample),
        epoch_of_next_sample(epochs_per_sample),
        epochs_per_negative_sample(epochs_per_sample.size()),
        epoch_of_next_negative_sample(epochs_per_sample.size()) {
    std::size_t n_edges = epochs_per_sample.size();
    float nsr = 1.0 / negative_sample_rate;
    for (std::size_t i = 0; i < n_edges; i++) {
      epochs_per_negative_sample[i] = epochs_per_sample[i] * nsr;
      epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
    }
  }

  auto epoch_begin(std::size_t epoch) -> void {
    this->epoch = epoch;
  }

  // are we due to sample the edge on the current epoch?
  auto is_sample_edge(std::size_t edge) const -> bool {
    return epoch_of_next_sample[edge] <= epoch;
  }

  auto get_num_neg_samples(std::size_t edge) const -> std::size_t {
    return static_cast<std::size_t>(
        (epoch - epoch_of_next_negative_sample[edge]) /
          epochs_per_negative_sample[edge]);
  }

  void next_sample(std::size_t edge, std::size_t num_neg_samples) {
    // set the next epoch when this edge will be sampled
    epoch_of_next_sample[edge] += epochs_per_sample[edge];
    epoch_of_next_negative_sample[edge] +=
        num_neg_samples * epochs_per_negative_sample[edge];
  }

private:
  // how often to sample each edge
  std::vector<float> epochs_per_sample;
  // the epoch when the edge should be sampled next
  std::vector<float> epoch_of_next_sample;
  std::vector<float> epochs_per_negative_sample;
  std::vector<float> epoch_of_next_negative_sample;
};


class NCVisSampler {
public:
  std::size_t epoch;
  NCVisSampler(const std::vector<float> &epochs_per_sample,
               const std::vector<std::size_t> &negative_plan)
    : epochs_per_sample(epochs_per_sample),
      epoch_of_next_sample(epochs_per_sample),
      negative_plan(negative_plan) {}

  auto epoch_begin(std::size_t epoch) -> void {
    this->epoch = epoch;
  }

  auto is_sample_edge(std::size_t edge) const -> bool {
    return epoch_of_next_sample[edge] <= epoch;
  }

  auto get_num_neg_samples(std::size_t edge) const -> std::size_t {
    return negative_plan[epoch];
  }

  void next_sample(std::size_t edge, std::size_t num_neg_samples) {
    epoch_of_next_sample[edge] += epochs_per_sample[edge];
  }

private:
  std::vector<float> epochs_per_sample;
  std::vector<float> epoch_of_next_sample;
  std::vector<std::size_t> negative_plan;
};

} // namespace uwot

#endif // UWOT_SAMPLER_H

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
  Sampler(const std::vector<float> &epochs_per_sample,
          float negative_sample_rate)

      :

        epochs_per_sample(epochs_per_sample),
        epoch_of_next_sample(epochs_per_sample),
        epochs_per_negative_sample(epochs_per_sample.size()),
        epoch_of_next_negative_sample(epochs_per_sample.size()) {
    std::size_t esz = epochs_per_sample.size();
    float nsr = 1.0 / negative_sample_rate;
    for (std::size_t i = 0; i < esz; i++) {
      epochs_per_negative_sample[i] = epochs_per_sample[i] * nsr;
      epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
    }
  }
  auto is_sample_edge(std::size_t i, std::size_t n) const -> bool {
    return epoch_of_next_sample[i] <= n;
  }
  auto get_num_neg_samples(std::size_t i, std::size_t n) const -> std::size_t {
    return static_cast<std::size_t>((n - epoch_of_next_negative_sample[i]) /
                                    epochs_per_negative_sample[i]);
  }

  void next_sample(std::size_t i, std::size_t num_neg_samples) {
    epoch_of_next_sample[i] += epochs_per_sample[i];
    epoch_of_next_negative_sample[i] +=
        num_neg_samples * epochs_per_negative_sample[i];
  }

private:
  std::vector<float> epochs_per_sample;
  std::vector<float> epoch_of_next_sample;
  std::vector<float> epochs_per_negative_sample;
  std::vector<float> epoch_of_next_negative_sample;
};

} // namespace uwot

#endif // UWOT_SAMPLER_H
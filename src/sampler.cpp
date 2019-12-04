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

#include "sampler.h"

Sampler::Sampler(const std::vector<double> &epochs_per_sample,
                 const double negative_sample_rate)
    :

      epochs_per_sample(epochs_per_sample),
      epoch_of_next_sample(epochs_per_sample),
      epochs_per_negative_sample(epochs_per_sample.size()),
      epoch_of_next_negative_sample(epochs_per_sample.size()) {
  const std::size_t esz = epochs_per_sample.size();
  const double nsr = 1.0 / negative_sample_rate;
  for (std::size_t i = 0; i < esz; i++) {
    epochs_per_negative_sample[i] = epochs_per_sample[i] * nsr;
    epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
  }
}

bool Sampler::is_sample_edge(const std::size_t i, const std::size_t n) const {
  return epoch_of_next_sample[i] <= n;
}

const std::size_t Sampler::get_num_neg_samples(const std::size_t i,
                                               const std::size_t n) const {
  return static_cast<std::size_t>((n - epoch_of_next_negative_sample[i]) /
                                  epochs_per_negative_sample[i]);
}

void Sampler::next_sample(const std::size_t i,
                          const std::size_t num_neg_samples) {
  epoch_of_next_sample[i] += epochs_per_sample[i];
  epoch_of_next_negative_sample[i] +=
      num_neg_samples * epochs_per_negative_sample[i];
}

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

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "sampler.h"

Sampler::Sampler(const arma::vec& epochs_per_sample, 
                 double negative_sample_rate) :

  _epochs_per_sample(epochs_per_sample.size()),
  epoch_of_next_sample(epochs_per_sample.size()),
  epochs_per_negative_sample(epochs_per_sample.size()),
  epoch_of_next_negative_sample(epochs_per_sample.size())
  {
    for (std::size_t i = 0; i < epochs_per_sample.size(); i++) {
      _epochs_per_sample[i] = epochs_per_sample[i];
      epoch_of_next_sample[i] = epochs_per_sample[i];
      epochs_per_negative_sample[i] = epochs_per_sample[i] / negative_sample_rate;
      epoch_of_next_negative_sample[i] = epochs_per_negative_sample[i];
    }
  }

bool Sampler::is_sample_edge(std::size_t i, std::size_t n) const {
  return epoch_of_next_sample[i] <= n;
}

unsigned int Sampler::get_num_neg_samples(std::size_t i, std::size_t n) const {
  auto n_neg_samples = static_cast<unsigned int>(
    (n - epoch_of_next_negative_sample[i]) / 
      epochs_per_negative_sample[i]);
  
  return n_neg_samples;
}

void Sampler::next_sample(std::size_t i, unsigned int num_neg_samples) {
  epoch_of_next_sample[i] += _epochs_per_sample[i];
  
  epoch_of_next_negative_sample[i] += 
    num_neg_samples * epochs_per_negative_sample[i];
}

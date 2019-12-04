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

#ifndef UWOT_SAMPLER_H
#define UWOT_SAMPLER_H

#include <vector>

// Weighted edge sampler
class Sampler {
public:
  Sampler(const std::vector<double> &epochs_per_sample,
          const double negative_sample_rate);

  bool is_sample_edge(const std::size_t i, const std::size_t n) const;
  const std::size_t get_num_neg_samples(const std::size_t i,
                                        const std::size_t n) const;
  void next_sample(const std::size_t i, const std::size_t num_neg_samples);

private:
  std::vector<double> epochs_per_sample;
  std::vector<double> epoch_of_next_sample;
  std::vector<double> epochs_per_negative_sample;
  std::vector<double> epoch_of_next_negative_sample;
};

#endif // UWOT_SAMPLER_H
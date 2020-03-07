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

#ifndef UWOT_RNG_H
#define UWOT_RNG_H

#include <limits>

// linked from dqrng
#include "convert_seed.h"
#include "pcg_random.hpp"

#include "uwot/tauprng.h"

// NOT THREAD SAFE
// based on code in the dqsample package
static uint64_t random64() {
  return R::runif(0, 1) * (std::numeric_limits<uint64_t>::max)();
}

// NOT THREAD SAFE
static uint32_t random32() {
  return R::runif(0, 1) * (std::numeric_limits<uint32_t>::max)();
}

struct tau_factory {
  uint64_t seed1;
  uint64_t seed2;
  tau_factory() : seed1(random64()), seed2(random64()) {}

  void reseed() {
    seed1 = random64();
    seed2 = random64();
  }

  uwot::tau_prng create(uint64_t seed) {
    return uwot::tau_prng(seed1, seed2, seed);
  }
};

struct pcg_prng {

  pcg32 gen;

  pcg_prng(uint64_t seed) { gen.seed(seed); }

  // return a value in (0, n]
  std::size_t operator()(std::size_t n) {
    std::size_t result = gen(n);
    return result;
  }
};

struct pcg_factory {
  uint32_t seed1;
  pcg_factory() : seed1(random32()) {}

  void reseed() { seed1 = random32(); }

  pcg_prng create(uint32_t seed) {
    uint32_t seeds[2] = {seed1, seed};
    return pcg_prng(dqrng::convert_seed<uint64_t>(seeds, 2));
  }
};

#endif // UWOT_RNG_H

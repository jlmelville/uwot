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

// Three-component combined Tausworthe "taus88" PRNG from L'Ecuyer 1996.

#ifndef UWOT_TAUPRNG_H
#define UWOT_TAUPRNG_H

#include "Rcpp.h"
#include <limits>
// linked from dqrng
#include "convert_seed.h"
#include "pcg_random.hpp"

// NOT THREAD SAFE
// based on code in the dqsample package
static uint64_t random64() {
  return R::runif(0, 1) * (std::numeric_limits<uint64_t>::max)();
}

// NOT THREAD SAFE
static uint32_t random32() {
  return R::runif(0, 1) * (std::numeric_limits<uint32_t>::max)();
}

struct tau_prng {
  uint64_t state0;
  uint64_t state1; // technically this needs to always be > 7
  uint64_t state2; // and this should be > 15

  static constexpr uint64_t MAGIC0 = static_cast<uint64_t>(4294967294);
  static constexpr uint64_t MAGIC1 = static_cast<uint64_t>(4294967288);
  static constexpr uint64_t MAGIC2 = static_cast<uint64_t>(4294967280);

  tau_prng(uint64_t state0, uint64_t state1, uint64_t state2)
      : state0(state0), state1(state1 > 7 ? state1 : 8),
        state2(state2 > 15 ? state2 : 16) {}

  int32_t operator()() {
    state0 = (((state0 & MAGIC0) << 12) & 0xffffffff) ^
             ((((state0 << 13) & 0xffffffff) ^ state0) >> 19);
    state1 = (((state1 & MAGIC1) << 4) & 0xffffffff) ^
             ((((state1 << 2) & 0xffffffff) ^ state1) >> 25);
    state2 = (((state2 & MAGIC2) << 17) & 0xffffffff) ^
             ((((state2 << 3) & 0xffffffff) ^ state2) >> 11);

    return state0 ^ state1 ^ state2;
  }

  // return a value in (0, n]
  std::size_t operator()(const std::size_t n) {
    std::size_t result = (*this)() % n;
    return result;
  }
};

struct tau_factory {
  uint64_t seed1;
  uint64_t seed2;
  tau_factory() : seed1(random64()), seed2(random64()) {}

  void reseed() {
    seed1 = random64();
    seed2 = random64();
  }

  tau_prng create(uint64_t seed) { return tau_prng(seed1, seed2, seed); }
};

struct pcg_prng {

  pcg32 gen;

  pcg_prng(uint64_t seed) { gen.seed(seed); }

  // return a value in (0, n]
  std::size_t operator()(const std::size_t n) {
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

#endif // UWOT_TAUPRNG_H

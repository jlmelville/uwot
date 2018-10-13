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

class tau_prng {
  long long state0;
  long long state1;
  long long state2;
public:
  tau_prng(long long state0, long long state1, long long state2)
    : state0(state0), state1(state1), state2(state2) {}

  int operator()() {
    state0 = (((state0 & 4294967294LL) << 12) & 0xffffffff) ^
      ((((state0 << 13) & 0xffffffff) ^ state0) >> 19);
    state1 = (((state1 & 4294967288LL) << 4) & 0xffffffff) ^
      ((((state1 << 2) & 0xffffffff) ^ state1) >> 25);
    state2 = (((state2 & 4294967280LL) << 17) & 0xffffffff) ^
      ((((state2 << 3) & 0xffffffff) ^ state2) >> 11);

    return state0 ^ state1 ^ state2;
  }
};

#endif // UWOT_TAUPRNG_H

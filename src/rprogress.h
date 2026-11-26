//  UWOT -- An R package for dimensionality reduction using UMAP
//
//  Copyright (C) 2021 James Melville
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
//  along with UWOT.  If not, see <http://www.gnu.org/licenses

#ifndef UWOT_RPROGRESS_H
#define UWOT_RPROGRESS_H

#include <interruptable_progress_monitor.hpp>
#include <simple_progress_bar.hpp>

struct RProgress {
  SimpleProgressBar bar;
  InterruptableProgressMonitor monitor;

  RProgress(std::size_t n_epochs, bool verbose)
      : monitor(n_epochs, verbose, bar) {}

  auto is_aborted() -> bool { return monitor.check_abort(); }

  void report() { monitor.increment(); }
};

#endif // UWOT_RPROGRESS_H

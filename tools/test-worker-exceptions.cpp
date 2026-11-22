#include "pforr.h"

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <string>

template <typename Run> bool check_failure(Run run, std::size_t threads) {
  std::atomic<std::size_t> completed{0};
  auto worker = [&](std::size_t begin, std::size_t end) {
    if (begin == 0) {
      throw std::logic_error("worker failed");
    }
    completed += end - begin;
  };

  try {
    run(worker, threads);
  } catch (const std::exception &error) {
    return std::string(error.what()) == "worker failed" &&
           completed == (threads == 1 ? 0 : 8) &&
           (threads != 1 || dynamic_cast<const std::logic_error *>(&error));
  }
  return false;
}

int main() {
  auto regular = [](auto &worker, std::size_t threads) {
    pforr::parallel_for(0, 16, worker, threads, 1);
  };
  auto indexed = [](auto &worker, std::size_t threads) {
    auto adapter = [&](std::size_t begin, std::size_t end, std::size_t) {
      worker(begin, end);
    };
    pforr::parallel_for_indexed(0, 16, adapter, threads, 1);
  };
  for (std::size_t threads : {1, 2}) {
    if (!check_failure(regular, threads) || !check_failure(indexed, threads)) {
      std::cerr << "Worker exception check failed with " << threads
                << " threads\n";
      return 1;
    }
  }
}

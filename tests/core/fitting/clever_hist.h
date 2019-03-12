#pragma once

#include <vector>
#include <cstddef>

struct CleverHist
{
  std::vector<double> bins;
  std::vector<double> counts;

  CleverHist() = default;

  static CleverHist make_linear(const std::vector<double>& events, size_t bin_count);
};
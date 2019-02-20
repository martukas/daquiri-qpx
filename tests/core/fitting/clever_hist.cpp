#include "clever_hist.h"

#include <algorithm>
#include <cmath>

CleverHist CleverHist::make_linear(const std::vector<double>& events, size_t bin_count)
{
  CleverHist ret;
  if (events.empty() || !bin_count)
    return ret;

  double min = *std::min_element(events.begin(), events.end());
  double max = *std::max_element(events.begin(), events.end());

  double range = max - min;
  double bin_width = (range) / static_cast<double>(bin_count);
  range += 0.01 * bin_width;
  bin_width = (range) / static_cast<double>(bin_count);

  ret.bins.resize(bin_count, 0.0);
  ret.counts.resize(bin_count, 0.0);

  for (size_t i=0; i < ret.bins.size(); ++i)
    ret.bins[i] = min + i * bin_width + 0.5 * bin_width;

  for (const auto& e : events)
    ret.counts[static_cast<size_t>((e - min)/ bin_width)]++;

  return ret;
}

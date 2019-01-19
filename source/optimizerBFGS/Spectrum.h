#pragma once

#include <optimizerBFGS/Calibration.h>
#include <vector>

namespace Hypermet
{

class CSpectrum
{
 public:
  std::vector<size_t> channels;
  Calibration calibration;
  double weight(size_t i) const;

  static double dead_time(double real_time, double live_time);
  static double rate(double live_time, double sum_counts);

 private:
  // \todo template type for val;
  size_t mystery_function(double val);
};

}

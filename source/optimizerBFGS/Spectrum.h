#pragma once

#include <optimizerBFGS/Calibration.h>
#include <vector>

namespace Hypermet
{

class CSpectrum
{
 public:
  std::vector<size_t> Channel;
  Calibration calibration;
  double Weight(size_t i) const;

  static double DeadTime(double TrueTime, double LiveTime);
  static double Rate(double LiveTime, double SumCounts);

 private:
  // \todo template type for val;
  size_t mystery_function(double val);
};

}

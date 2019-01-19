#pragma once

#include <optimizerBFGS/Calibration.h>
#include <vector>

namespace Hypermet
{

class CSpectrum
{
 public:
  std::vector<size_t> Channel;
  CCalibration Calibration;
  double Weight(size_t i) const;

  static double DeadTime(double TrueTime, double LiveTime);
  static double Rate(double LiveTime, double SumCounts);

 private:
  // \todo template type for Val;
  size_t mystery_function(double Val);
};

}

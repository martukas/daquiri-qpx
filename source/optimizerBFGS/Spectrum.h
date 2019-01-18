#pragma once

#include <optimizerBFGS/Calibration.h>
#include <vector>

class CSpectrum
{
 public:
  std::vector<size_t> Channel;
  CCalibration Calibration;
  double Weight(size_t i) const;

  // template type for Val;
  int8_t Sign(double Val);

  double DeadTime(double TrueTime, double LiveTime);
  double Rate(double LiveTime, double SumCounts);

 private:
  // template type for Val;
  size_t mystery_function(double Val);
};


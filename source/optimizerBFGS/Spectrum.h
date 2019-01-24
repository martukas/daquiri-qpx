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
  double weight_true(size_t i) const;
  double weight_phillips_marlow(size_t i) const;
  double weight_revay_student(size_t i) const;

 private:
};


// \todo these go somewhere else

double dead_time(double real_time, double live_time);
double rate(double live_time, double sum_counts);
size_t mystery_function(double val);

}

#pragma once

#include <optimizerBFGS/Value.h>
#include <optimizerBFGS/Calibration.h>

namespace Hypermet
{

class Peak
{
 private:
  int32_t FEP_status_{1};
 public:
  ValueDefault position;
  ValueGam amplitude;

  int32_t step_type() const;
  double peak_position() const;
  double peak_position_unc() const;
  double peak_energy(const Calibration& cal) const;
  double peak_energy_unc(const Calibration& cal) const;
  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Peak& other) const;
};

}

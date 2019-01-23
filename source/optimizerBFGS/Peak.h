#pragma once

#include <optimizerBFGS/Step.h>
#include <optimizerBFGS/Tail.h>
#include <optimizerBFGS/Calibration.h>

namespace Hypermet
{

class Peak
{
 public:
  struct Components
  {
    double gaussian{0};
    double short_tail{0};
    double right_tail{0};
    double long_tail{0};
    double step{0};
  };

  // These are unique to peak
  Value position;
  ValueGam amplitude;

  // By default these are not unique to peak
  bool width_override{false};
  Value width_;

  // skews (part of peak)
  Tail short_tail {Side::left};
  Tail right_tail {Side::right};

  // step & tail (background)
  Tail long_tail {Side::left};
  Step step;

  double peak_position() const;
  double peak_position_unc() const;
  double peak_energy(const Calibration& cal) const;
  double peak_energy_unc(const Calibration& cal) const;
  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Peak& other) const;

  double area() const;
  double area_uncert(double chisq_norm) const;

  Components eval(double chan) const;
  Components eval_grad(double chan, std::vector<double>& grads) const;

  PrecalcVals precalc_vals(double chan) const;

  void update_indices(int32_t& i);

  void put(std::vector<double>& fit) const;
  void get(const std::vector<double>& fit);
  void get_uncerts(const std::vector<double>& diagonals, double chisq_norm);
};

}

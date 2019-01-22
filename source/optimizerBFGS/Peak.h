#pragma once

#include <optimizerBFGS/Value.h>
#include <optimizerBFGS/Calibration.h>

namespace Hypermet
{

struct PrecalcVals
{
  double width;
  double ampl;
  double half_ampl;
  double spread;
};


struct Tail
{
  bool override{false};
  bool enabled{true};
  Value amplitude, slope;

  static double eval_with(const PrecalcVals& pre, double ampl, double slp);
  static double eval_flipped_with(const PrecalcVals& pre, double ampl, double slp);

  double eval(const PrecalcVals& pre) const;
  double eval_flipped(const PrecalcVals& pre) const;

  double eval_grad(const PrecalcVals& pre,
                   std::vector<double>& grads,
                   size_t i_width, size_t i_pos, size_t i_amp);

  double eval_flipped_grad(const PrecalcVals& pre,
                            std::vector<double>& grads,
                            size_t i_width, size_t i_pos, size_t i_amp);
};

class Peak
{
 private:
  int32_t FEP_status_{1};

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

  // skews
  Tail short_tail;
  Tail right_tail;

  // step & tail
  Tail long_tail;

  bool step_override{false};
  bool step_enabled_{true};
  Value step_amplitude_;

  int32_t step_type() const;
  double peak_position() const;
  double peak_position_unc() const;
  double peak_energy(const Calibration& cal) const;
  double peak_energy_unc(const Calibration& cal) const;
  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Peak& other) const;

  Components eval(double chan) const;
  Components eval_grad(double chan, std::vector<double>& grads, size_t offset);

  PrecalcVals precalc_vals(double chan) const;

  static double eval_skew(double ampl, double spread, double slope);
};

}

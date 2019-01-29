#pragma once

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/hypermet/Tail.h>
#include <core/fitting/hypermet/Calibration.h>
#include <core/calibration/calibration.h>
#include <core/fitting/uncertain.h>

namespace DAQuiri
{

class Hypermet
{
 public:
  struct Components
  {
    double gaussian{0};
    double short_tail{0};
    double right_tail{0};
    double long_tail{0};
    double step{0};

    double peak_skews() const;
    double step_tail() const;
    double all() const;
  };

  Hypermet();
  void apply_defaults(const Hypermet& other);
  void force_defaults(const Hypermet& other);

  Hypermet gaussian_only() const;
  bool is_gaussian_only() const;

  bool sanity_check(double min_x, double max_x) const;

  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Hypermet& other) const;

  void update_indices(int32_t& i);
  void put(std::vector<double>& fit) const;
  void get(const std::vector<double>& fit);
  void get_uncerts(const std::vector<double>& diagonals, double chisq_norm);

  UncertainDouble peak_position() const;
  UncertainDouble peak_energy(const HCalibration& cal) const;
  UncertainDouble peak_energy(const Calibration& cal) const;
  UncertainDouble area() const;
  UncertainDouble peak_area_eff(const HCalibration& cal) const;

  PrecalcVals precalc_vals(double chan) const;
  PrecalcVals precalc_vals_at(double chan, const std::vector<double>& fit) const;

  Components eval(double chan) const;
  Components eval_at(double chan, const std::vector<double>& fit) const;
  Components eval_grad(double chan, std::vector<double>& grads) const;
  Components eval_grad_at(double chan, const std::vector<double>& fit,
      std::vector<double>& grads) const;

  std::string to_string() const;

  // These are unique to peak
  Value position;
  ValueGam amplitude;

  // By default these are not unique to peak
  bool width_override{false};
  Value width_;

  // \todo why skew naming different?
  // skews (part of peak)
  Tail short_tail {Side::left};
  Tail right_tail {Side::right};

  // step & tail (background)
  Tail long_tail {Side::left};
  Step step;

  double chi_sq_norm {0.0};
};

void to_json(nlohmann::json& j, const Hypermet& s);
void from_json(const nlohmann::json& j, Hypermet& s);

}

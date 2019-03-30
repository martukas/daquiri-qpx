#pragma once

#include <core/fitting/hypermet/Step.h>
#include <core/fitting/hypermet/Tail.h>
#include <core/fitting/hypermet/Calibration.h>
#include <core/calibration/calibration.h>
#include <core/fitting/sum4/sum4.h>

namespace DAQuiri
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

    double peak_skews() const;
    double step_tail() const;
    double all() const;
  };

  Peak();
  // \todo set position bounds
  void apply_defaults(const Peak& other);
  void force_defaults(const Peak& other);

  Peak gaussian_only() const;
  bool is_gaussian_only() const;

  bool sanity_check(double min_x, double max_x) const;
  bool sane() const;

  bool full_energy_peak() const;
  void full_energy_peak(bool flag);
  bool operator<(const Peak& other) const;

  void reset_indices();
  void update_indices(int32_t& i);
  void put(Eigen::VectorXd& fit) const;
  void get(const Eigen::VectorXd& fit);
  void get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm);

  double id() const;

  UncertainDouble peak_position() const;
  UncertainDouble peak_energy(const HCalibration& cal) const;
  UncertainDouble peak_energy(const Calibration& cal) const;

  UncertainDouble fwhm() const;
  UncertainDouble fwhm_energy(const HCalibration& cal) const;
  UncertainDouble fwhm_energy(const Calibration& cal) const;

  UncertainDouble area() const;
  // \todo use duration type with uncertainty
  UncertainDouble rate(double live_time) const;
  // \todo get rate rather than area
  UncertainDouble peak_area_eff(const HCalibration& cal) const;

  PrecalcVals precalc_vals(double chan) const;
  PrecalcVals precalc_vals_at(double chan, const Eigen::VectorXd& fit) const;

  Components eval(double chan) const;
  Components eval_at(double chan, const Eigen::VectorXd& fit) const;
  Components eval_grad(double chan, Eigen::VectorXd& grads) const;
  Components eval_grad_at(double chan, const Eigen::VectorXd& fit, Eigen::VectorXd& grads) const;

  std::string to_string(std::string prepend = "") const;

  // These are unique to peak
  SineBoundedValue position;
  PositiveValue amplitude; // was ValueGam

  // By default these are not unique to peak
  bool width_override{false};
  SineBoundedValue width;

  // \todo why skew naming different?
  // skews (part of peak)
  Tail short_tail {Side::left};
  Tail right_tail {Side::right};

  // step & tail (background)
  Tail long_tail {Side::left};
  Step step;

  double chi_sq_norm {0.0};

  SUM4 sum4;

  // \todo these might have to be somewhere else
  UncertainDouble true_energy;
  UncertainDouble theoretical_intensity {0.0, 0.0};
  UncertainDouble relative_efficiency {0.0, 0.0};
};

void to_json(nlohmann::json& j, const Peak& s);
void from_json(const nlohmann::json& j, Peak& s);

}

#pragma once

#include <core/fitting/gaussian.h>
#include <core/fitting/fit_settings.h>

namespace DAQuiri
{

class Hypermet
{
 public:
  Hypermet() :
      Hypermet(Gaussian(), FitSettings()) {}

  Hypermet(Gaussian gauss, FitSettings settings);

  const Parameter& center() const { return center_; }
  const Parameter& height() const { return height_; }
  const Parameter& width() const { return width_; }
  const Parameter& Lskew_amplitude() const { return Lskew_amplitude_; }
  const Parameter& Lskew_slope() const { return Lskew_slope_; }
  const Parameter& Rskew_amplitude() const { return Rskew_amplitude_; }
  const Parameter& Rskew_slope() const { return Rskew_slope_; }
  const Parameter& tail_amplitude() const { return tail_amplitude_; }
  const Parameter& tail_slope() const { return tail_slope_; }
  const Parameter& step_amplitude() const { return step_amplitude_; }
  double chi2() const { return chi2_; }

  bool user_modified() const { return user_modified_; }

  // \todo overload with uncertain type for these
  void set_center(const Parameter& ncenter);
  void set_height(const Parameter& nheight);
  void set_width(const Parameter& nwidth);
  void set_Lskew_amplitude(const Parameter& nLskew_amplitude);
  void set_Lskew_slope(const Parameter& nLskew_slope);
  void set_Rskew_amplitude(const Parameter& nRskew_amplitude);
  void set_Rskew_slope(const Parameter& nRskew_slope);
  void set_tail_amplitude(const Parameter& ntail_amplitude);
  void set_tail_slope(const Parameter& ntail_slope);
  void set_step_amplitude(const Parameter& nstep_amplitude);

  void set_chi2(double);

  void constrain_center(double min, double max);
  void constrain_height(double min, double max);
  void constrain_width(double min, double max);

  std::string to_string() const;
  // \todo use uncertain type
  double eval_peak(double) const;
  double eval_step_tail(double) const;
  std::vector<double> peak(std::vector<double> x) const;
  std::vector<double> step_tail(std::vector<double> x) const;
  // \todo use uncertain type
  double area() const;
  bool gaussian_only() const;
  Gaussian gaussian() const;

  friend void to_json(nlohmann::json& j, const Hypermet& s);
  friend void from_json(const nlohmann::json& j, Hypermet& s);

 private:
  Parameter height_, center_, width_,
      Lskew_amplitude_, Lskew_slope_,
      Rskew_amplitude_, Rskew_slope_,
      tail_amplitude_, tail_slope_,
      step_amplitude_;

  bool step_enabled_{true};
  bool tail_enabled_{true};
  bool Lskew_enabled_{true};
  bool Rskew_enabled_{true};

  double chi2_{0};
  bool user_modified_{false};
};

}

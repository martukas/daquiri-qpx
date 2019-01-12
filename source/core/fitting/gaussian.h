#pragma once

#include <core/calibration/polynomial.h>

namespace DAQuiri
{

class Gaussian
{
 public:
  Gaussian() = default;

  std::string to_string() const;

  // \todo use uncertain
  double evaluate(double x);
  std::vector<double> evaluate_array(std::vector<double> x);

  // \todo use uncertain
  double area() const;

  const Parameter& center() const { return center_; }
  const Parameter& height() const { return height_; }
  const Parameter& hwhm() const { return hwhm_; }

  // \todo overload with uncertain type for these 3
  void set_center(const Parameter& ncenter);
  void set_height(const Parameter& nheight);
  void set_hwhm(const Parameter& nwidth);

  void constrain_center(double min, double max);
  void constrain_height(double min, double max);
  void constrain_hwhm(double min, double max);

  void set_chi2(double);
  // \todo get chi2

 protected:
  Parameter height_{0};
  Parameter center_{0};
  Parameter hwhm_{0};
  double chi2_{0};
};

}
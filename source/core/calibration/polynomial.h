#pragma once

#include <core/calibration/coef_function.h>
#include <vector>
#include <core/fitting/hypermet/Value.h>

namespace DAQuiri
{

class Polynomial : public CoefFunction
{
 public:
  using CoefFunction::CoefFunction;

  Polynomial(const std::vector<double>& coeffs, double uncert = 0.0);

  bool valid() const override;
  std::string type() const override { return "Polynomial"; }
  Polynomial* clone() const override { return new Polynomial(*this); }
  bool is_equal(CoefFunction* other) const override;
  double operator() (double x) const override;
  double derivative(double) const override;

  std::string debug() const override;
  std::string to_UTF8(int precision, bool with_rsq) const override;
  std::string to_markup(int precision, bool with_rsq) const override;

 protected:
  std::vector<UnboundedValue> coeffs_;
};

}
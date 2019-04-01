#pragma once

#include <core/calibration/coef_function.h>
#include <core/fitting/data_model/value.h>
#include <vector>

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
  double eval(double x) const override;
  double d_dx(double x) const override;

  void update_indices() override;
  Eigen::VectorXd variables() const override;
  double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override;
  void save_fit(const DAQuiri::FitResult& result) override;

  std::string to_string(std::string prepend = "") const override;
  std::string to_UTF8(int precision, bool with_rsq) const override;
  std::string to_markup(int precision, bool with_rsq) const override;

 protected:
  std::vector<UnboundedValue> coeffs_;
};

}
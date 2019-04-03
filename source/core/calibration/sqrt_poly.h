#pragma once

#include <core/calibration/calib_function.h>
#include <core/fitting/parameter/unbounded_param.h>
#include <vector>

namespace DAQuiri
{

class SqrtPoly : public CalibFunction
{
 public:
  using CalibFunction::CalibFunction;

  SqrtPoly(const std::vector<double>& coeffs);

  bool valid() const override;
  std::string type() const override { return "SqrtPoly"; }
  SqrtPoly* clone() const override { return new SqrtPoly(*this); }
  bool is_equal(CalibFunction* other) const override;
  double eval(double x) const override;
  double d_dx(double x) const override;

  void update_indices() override;
  Eigen::VectorXd variables() const override;
  double eval_grad_at(double x, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override;
  void save_fit(const DAQuiri::FitResult& result) override;

  std::string to_string(std::string prepend = "") const override;
  std::string to_UTF8(int precision) const override;
  std::string to_markup(int precision) const override;

  nlohmann::json to_json() const override {};
  void from_json(const nlohmann::json&) override {};

 protected:
  std::vector<UnboundedParam> coeffs_;
};

}
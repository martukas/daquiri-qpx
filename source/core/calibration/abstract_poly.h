#pragma once

#include <core/calibration/function.h>
#include <core/fitting/parameter/unbounded_param.h>
#include <vector>

namespace DAQuiri
{

class AbstractPolynomial : public CalibFunction
{
 public:
  using CalibFunction::CalibFunction;

  std::vector<UnboundedParam> coefficients;

  AbstractPolynomial() = default;
  AbstractPolynomial(const std::vector<double>& coeffs);

  bool valid() const override;
  bool is_equal(CalibFunction* other) const override;

  void update_indices() override;
  Eigen::VectorXd variables() const override;
  void save_fit(const DAQuiri::FitResult& result) override;

  nlohmann::json to_json() const override;
  void from_json(const nlohmann::json&) override;

  std::string to_string(std::string prepend = "") const override;
};

}
#pragma once

#include <core/calibration/polynomial.h>

namespace DAQuiri
{

class SqrtPoly : public Polynomial
{
 public:
  using Polynomial::Polynomial;

  std::string type() const override { return "SqrtPoly"; }
  SqrtPoly* clone() const override { return new SqrtPoly(*this); }

  double eval(double x) const override;
  double eval_at(double x, const Eigen::VectorXd& fit) const override;
  double eval_grad_at(double x, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const override;
  double d_dx(double x) const override;

  std::string to_string(std::string prepend = "") const override;
  std::string to_UTF8(int precision) const override;
  std::string to_markup(int precision) const override;
};

}
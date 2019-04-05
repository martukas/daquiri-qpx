#pragma once

#include <core/calibration/abstract_poly.h>

namespace DAQuiri
{

class Polynomial : public AbstractPolynomial
{
 public:
  using AbstractPolynomial::AbstractPolynomial;

  std::string type() const override { return "Polynomial"; }
  Polynomial* clone() const override { return new Polynomial(*this); }

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
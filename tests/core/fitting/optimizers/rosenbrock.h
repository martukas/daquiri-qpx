#pragma once

#include <core/fitting/optimizers/fittable_function.h>

class Rosenbrock : public DAQuiri::FittableFunction
{
 private:
  size_t n_;

 public:
  Eigen::VectorXd vals_;

  Rosenbrock(size_t n)
      : n_(n)
  {
    vals_ = Eigen::VectorXd::Zero(n_);
  }

  Eigen::VectorXd variables() const override
  {
    return vals_;
  }

  double chi_sq(const Eigen::VectorXd& fit) const override
  {
    Eigen::VectorXd throw_away;
    return chi_sq_gradient(fit, throw_away);
  }

  double chi_sq_gradient(const Eigen::VectorXd& fit,
                         Eigen::VectorXd& gradients) const override
  {
    double fx = 0.0;
    gradients.setConstant(n_, 0.0);
    for(size_t i = 0; i < n_; i += 2)
    {
      double t1 = 1.0 - fit[i];
      double t2 = 10 * (fit[i + 1] - fit[i] * fit[i]);
      gradients[i + 1] = 20 * t2;
      gradients[i]     = -2.0 * (fit[i] * gradients[i + 1] + t1);
      fx += t1 * t1 + t2 * t2;
    }
    return fx;
  }
};

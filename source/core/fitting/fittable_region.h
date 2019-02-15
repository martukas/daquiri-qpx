#pragma once

#include <core/fitting/weighted_data.h>
#include <core/fitting/fittable_function.h>

namespace DAQuiri
{

class FittableRegion : public Fittable
{
 public:
  int32_t var_count{0};
  WeightedData data;

  virtual double eval(double chan) const = 0;

  virtual double eval_at(double chan, const Eigen::VectorXd& fit) const = 0;

  virtual double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                              Eigen::VectorXd& grads) const = 0;

  double chi_sq(const Eigen::VectorXd& fit) const override;

  double chi_sq_gradient(const Eigen::VectorXd& fit,
                         Eigen::VectorXd& gradients) const override;

  double degrees_of_freedom() const;
};

}
#pragma once

#include <core/fitting/optimizers/fittable_function.h>
#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

class FittableRegion : public FittableFunction
{
 public:
  int32_t variable_count{0};
  WeightedData data;

  FittableRegion() = default;
  virtual ~FittableRegion() = default;

  virtual double eval(double chan) const = 0;

  virtual double eval_at(double chan, const Eigen::VectorXd& fit) const = 0;

  virtual double eval_grad_at(double chan, const Eigen::VectorXd& fit,
                              Eigen::VectorXd& grads) const = 0;

  double chi_sq() const;

  //Calculates the Chi-square over a region
  double chi_sq(const Eigen::VectorXd& fit) const override;

  //Calculates the Chi-square and its gradient
  double chi_sq_gradient(const Eigen::VectorXd& fit,
                         Eigen::VectorXd& gradients) const override;

  double degrees_of_freedom() const;
};

}
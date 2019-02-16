#include <core/fitting/fittable_region.h>
#include <core/util/more_math.h>

namespace DAQuiri
{

double FittableRegion::degrees_of_freedom() const
{
  double dof = static_cast<double>(data.data.size()) - static_cast<double>(variable_count);
  if (dof < 0.0)
    return 0.0;
  return dof;
}

double FittableRegion::chi_sq() const
{
  double chi_squared {0.0};
  for (const auto& p : data.data)
    chi_squared += square((p.count - this->eval(p.channel)) / p.weight_phillips_marlow);
  return chi_squared / degrees_of_freedom();
}

double FittableRegion::chi_sq(const Eigen::VectorXd& fit) const
{
  double chi_squared {0.0};
  for (const auto& p : data.data)
    chi_squared += square((p.count - this->eval_at(p.channel, fit)) / p.weight_phillips_marlow);
  return chi_squared / degrees_of_freedom();
}

double FittableRegion::chi_sq_gradient(const Eigen::VectorXd& fit,
                                       Eigen::VectorXd& gradients) const
{
  gradients.setConstant(fit.size(), 0.0);

  double chi_squared{0.0};
  Eigen::VectorXd channel_gradients;
  for (const auto& p : data.data)
  {
    channel_gradients.setConstant(fit.size(), 0.0);

    double fit_value = this->eval_grad_at(p.channel, fit, channel_gradients);
    chi_squared += square((p.count - fit_value) / p.weight_phillips_marlow);

    double grad_norm = -2.0 * (p.count - fit_value) / square(p.weight_phillips_marlow);
    for (size_t var = 0; var < fit.size(); ++var)
      gradients[var] += channel_gradients[var] * grad_norm;
  }
  return chi_squared / degrees_of_freedom();
}

}
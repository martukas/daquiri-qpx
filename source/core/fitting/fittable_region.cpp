#include <core/fitting/fittable_region.h>
#include <core/util/more_math.h>
#include <range/v3/all.hpp>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

double FittableRegion::degrees_of_freedom() const
{
  double dof = static_cast<double>(data.chan.size()) - static_cast<double>(variable_count);
  if (dof < 0.)
    return 0.;
  return dof;
}

double FittableRegion::chi_sq() const
{
  double chi_squared {0.};
  for (const auto& p : ranges::view::zip(data.chan, data.count, data.count_weight))
    chi_squared += square((std::get<1>(p) - this->eval(std::get<0>(p))) / std::get<2>(p));
  return chi_squared;
}

double FittableRegion::chi_sq(const Eigen::VectorXd& fit) const
{
  double chi_squared {0.};
  for (const auto& p : ranges::view::zip(data.chan, data.count, data.count_weight))
    chi_squared += square((std::get<1>(p) - this->eval_at(std::get<0>(p), fit))
        / std::get<2>(p));
  return chi_squared;
}

double FittableRegion::chi_sq_gradient(const Eigen::VectorXd& fit,
                                       Eigen::VectorXd& gradients) const
{
  gradients.setConstant(fit.size(), 0.);

  double chi_squared{0.};
  Eigen::VectorXd channel_gradients;
  for (const auto& p : ranges::view::zip(data.chan, data.count, data.count_weight))
  {
    channel_gradients.setConstant(fit.size(), 0.);

    double fit_value = this->eval_grad_at(std::get<0>(p), fit, channel_gradients);

    chi_squared += square((std::get<1>(p) - fit_value) / std::get<2>(p));

    double grad_norm = -2. * (std::get<1>(p) - fit_value) / square(std::get<2>(p));

    for (size_t var = 0; var < (size_t)fit.size(); ++var)
      gradients[var] += channel_gradients[var] * grad_norm;
  }
  return chi_squared;
}

}
#include <core/fitting/optimizers/BFGS.h>
#include <core/fitting/optimizers/brent.h>
#include <core/util/more_math.h>

#include <core/util/logger.h>

namespace DAQuiri
{

FitResult BudapestOptimizer::minimize(FittableFunction* fittable)
{
  FitResult ret;

  Brent brent;

  ret.variables = fittable->variables();
  auto var_count = static_cast<size_t>(ret.variables.size());
  Eigen::VectorXd
      search_direction(var_count),
      gradients(var_count),
      prev_val(var_count),
      Adg(var_count);

  double f = fittable->chi_sq_gradient(ret.variables, gradients);

  ret.inv_hessian.resize(var_count, var_count);
  ret.inv_hessian.setIdentity();

  for (; ret.iterations < maximum_iterations; ++ret.iterations)
  {
    double lambda = brent.line_search(fittable, ret.variables, search_direction);

    if (verbosity)
      DBG("lambda={}", lambda);

    ret.variables += lambda * search_direction;

    double fmin = fittable->chi_sq(ret.variables);

    //if (std::abs(f - fmin) < 0.000001) { done = true; }
    ret.converged = (2 * std::abs(fmin - f)) <= (tolerance * (std::abs(fmin) + std::abs(f) + eps));
    f = fmin;

    prev_val = gradients;

    fmin = fittable->chi_sq_gradient(ret.variables, gradients);
    if (verbosity)
      INFO("Fitting... Iteration = {}, Chisq = {}", ret.iterations, fmin);

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = gradients[i] - prev_val[i];

    for (size_t i = 0; i < var_count; ++i)
    {
      auto& Adgi = Adg[i];
      Adgi = 0.0;
      for (size_t j = 0; j < var_count; ++j)
        Adgi += ret.inv_hessian.coeff(i, j) * prev_val[j];
    }

    double s{0.0};
    double s1{0.0};
    for (size_t i = 0; i < var_count; ++i)
    {
      s1 += prev_val[i] * search_direction[i];
      s += prev_val[i] * Adg[i];
    }

    if (s1 != 0.0)
      s1 = 1.0 / s1;

    double s2{0.0};
    if (s != 0)
      s2 = 1.0 / s;

    for (size_t i = 0; i < var_count; ++i)
      prev_val[i] = s1 * search_direction[i] - s2 * Adg[i];

    for (size_t i = 0; i < var_count; ++i)
      for (size_t j = 0; j < var_count; ++j)
        ret.inv_hessian.coeffRef(i, j) +=
            s1 * square(search_direction[i])
                - s2 * square(Adg[i])
                + s * square(prev_val[i]);

    for (size_t i = 0; i < var_count; ++i)
    {
      search_direction[i] = 0;
      for (size_t j = 0; j < var_count; ++j)
        search_direction[i] -= ret.inv_hessian.coeff(i, j) * gradients[j];
    }

    if (ret.converged || cancel.load())
      break;
  }

  return ret;
}

}

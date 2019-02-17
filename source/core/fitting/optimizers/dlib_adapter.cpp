#include <core/fitting/optimizers/dlib_adapter.h>
#include <dlib/matrix/matrix_mat.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

FitResult DLibOptimizer::minimize(FittableFunction* fittable)
{
  function_ = fittable;

  variables_ = function_->variables();
  gradients_.setConstant(variables_.size(), 0.0);

  fitter_vector variables = dlib::mat(variables_);
  const auto& fe = std::bind(&DLibOptimizer::eval, this, std::placeholders::_1);
  const auto& fd = std::bind(&DLibOptimizer::derivative, this, std::placeholders::_1);

  FitResult ret;

  auto stop_strategy = dlib::objective_delta_stop_strategy(tolerance, maximum_iterations);
  if (verbose)
    stop_strategy = stop_strategy.be_verbose();

  ret.value = dlib::find_min(dlib::bfgs_search_strategy(),
                             stop_strategy, fe, fd, variables, minimum_value);

  ret.variables.setConstant(variables.size(), 0.0);
  for (long i = 0; i < variables.size(); ++i)
    ret.variables[i] = variables(i);

  return ret;
}

double DLibOptimizer::eval(const DLibOptimizer::fitter_vector& vars)
{
  for (long i = 0; i < vars.size(); ++i)
    variables_[i] = vars(i);
  return function_->chi_sq(variables_);
}

DLibOptimizer::fitter_vector DLibOptimizer::derivative(const DLibOptimizer::fitter_vector& vars)
{
  for (long i = 0; i < vars.size(); ++i)
    variables_[i] = vars(i);
  gradients_.setConstant(variables_.size(), 0.0);
  function_->chi_sq_gradient(variables_, gradients_);
  return dlib::mat(gradients_);
}

}

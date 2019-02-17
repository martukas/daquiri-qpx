#include <core/fitting/optimizers/gsl_adapter.h>

#include <memory>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

FitResult GSLOptimizer::minimize(FittableFunction* fittable)
{
  function_ = fittable;
  variables_ = function_->variables();
  gradients_.setConstant(variables_.size(), 0.0);

  gsl_multimin_function_fdf my_func;

  my_func.n = variables_.size();
  my_func.f = &GSLOptimizer::my_f;
  my_func.df = &GSLOptimizer::my_df;
  my_func.fdf = &GSLOptimizer::my_fdf;
  my_func.params = reinterpret_cast<void*>(this);

  std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>
      x(gsl_vector_alloc(variables_.size()), &gsl_vector_free);
  copy_to_gsl_vector(x.get(), variables_);

  std::unique_ptr<gsl_multimin_fdfminimizer, decltype(&gsl_multimin_fdfminimizer_free)>
      minimizer(gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, variables_.size()),
                &gsl_multimin_fdfminimizer_free);

  gsl_multimin_fdfminimizer_set(minimizer.get(), &my_func, x.get(), first_step_size, tolerance);

  FitResult ret;
  int status;
  do
  {
    ret.iterations++;
    status = gsl_multimin_fdfminimizer_iterate(minimizer.get());

    if (status)
      break;

    status = gsl_multimin_test_gradient(minimizer->gradient, gradient_tolerance);

    ret.variables = variables_;
    ret.value = minimizer->f;

    if (status == GSL_SUCCESS)
      ret.converged = true;

    if (verbose)
      INFO("{}", ret.to_string());
  }
  while (status == GSL_CONTINUE && ret.iterations < maximum_iterations);

  return ret;
}

double GSLOptimizer::my_f(const gsl_vector* v)
{
  copy_from_gsl_vector(variables_, v);
  return function_->chi_sq(variables_);
}

void GSLOptimizer::my_df(const gsl_vector* v, gsl_vector* df)
{
  copy_from_gsl_vector(variables_, v);
  gradients_.setConstant(variables_.size(), 0.0);
  function_->chi_sq_gradient(variables_, gradients_);
  copy_to_gsl_vector(df, gradients_);
}

void GSLOptimizer::my_fdf(double* ret, const gsl_vector* x, gsl_vector* df)
{
  copy_from_gsl_vector(variables_, x);
  gradients_.setConstant(variables_.size(), 0.0);
  *ret = function_->chi_sq_gradient(variables_, gradients_);
  copy_to_gsl_vector(df, gradients_);
}

double GSLOptimizer::my_f(const gsl_vector* v, void* params)
{
  auto this_adapter = reinterpret_cast<GSLOptimizer*>(params);
  return this_adapter->my_f(v);
}

void GSLOptimizer::my_df(const gsl_vector* v, void* params, gsl_vector* df)
{
  auto this_adapter = reinterpret_cast<GSLOptimizer*>(params);
  this_adapter->my_df(v, df);
}

void GSLOptimizer::my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
  auto this_adapter = reinterpret_cast<GSLOptimizer*>(params);
  this_adapter->my_fdf(f, x, df);
}

void GSLOptimizer::copy_to_gsl_vector(gsl_vector* gslv, const Eigen::VectorXd& eigenv)
{
  if (gslv->size != eigenv.size())
    throw std::runtime_error(
       fmt::format("copy_to_gsl_vector: vector size simatch {}!={}", gslv->size, eigenv.size()));
  for (size_t i = 0; i < eigenv.size(); ++i)
    gsl_vector_set(gslv, i, eigenv[i]);
}

void GSLOptimizer::copy_from_gsl_vector(Eigen::VectorXd& eigenv, const gsl_vector* gslv)
{
  if (gslv->size != eigenv.size())
    throw std::runtime_error(
        fmt::format("copy_from_gsl_vector: vector size simatch {}!={}", eigenv.size(), gslv->size));
  for (size_t i = 0; i < eigenv.size(); ++i)
    eigenv[i] = gsl_vector_get(gslv, i);
}

}

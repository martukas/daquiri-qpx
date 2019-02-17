#include <core/fitting/optimizers/gsl_adapter.h>

#include <memory>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

FitResult GSLAdapter::minimize(FittableFunction* fittable, double tolf)
{
  function_ = fittable;
  vars = function_->variables();
  grads.setConstant(vars.size(), 0.0);

  gsl_multimin_function_fdf my_func;

  my_func.n = vars.size();
  my_func.f = &GSLAdapter::my_f;
  my_func.df = &GSLAdapter::my_df;
  my_func.fdf = &GSLAdapter::my_fdf;
  my_func.params = reinterpret_cast<void*>(this);

  std::unique_ptr<gsl_vector, decltype(&gsl_vector_free)>
      x(gsl_vector_alloc(vars.size()), &gsl_vector_free);
  copy_to_gsl_vector(x.get(), vars);

  std::unique_ptr<gsl_multimin_fdfminimizer, decltype(&gsl_multimin_fdfminimizer_free)>
      minimizer(gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, vars.size()),
                &gsl_multimin_fdfminimizer_free);

  gsl_multimin_fdfminimizer_set(minimizer.get(), &my_func, x.get(), 0.01, 1e-4);

  FitResult ret;
  int status;
  do
  {
    ret.iterations++;
    status = gsl_multimin_fdfminimizer_iterate(minimizer.get());

    if (status)
      break;

    status = gsl_multimin_test_gradient(minimizer->gradient, 1e-3);

    if (status == GSL_SUCCESS)
      ret.converged = true;

    ret.variables = vars;
//    std::stringstream ss;
//    ss << ret.variables.transpose();
//    INFO("iter={}  vars={}  f={}", ret.iterations, ss.str(), minimizer->f);
  }
  while (status == GSL_CONTINUE && ret.iterations < 3000);

  ret.variables = vars;

  return ret;
}

double GSLAdapter::my_f(const gsl_vector* v)
{
  copy_from_gsl_vector(vars, v);
  return function_->chi_sq(vars);
}

void GSLAdapter::my_df(const gsl_vector* v, gsl_vector* df)
{
  copy_from_gsl_vector(vars, v);
  grads.setConstant(vars.size(), 0.0);
  function_->chi_sq_gradient(vars, grads);
  copy_to_gsl_vector(df, grads);
}

void GSLAdapter::my_fdf(double* ret, const gsl_vector* x, gsl_vector* df)
{
  copy_from_gsl_vector(vars, x);
  grads.setConstant(vars.size(), 0.0);
  *ret = function_->chi_sq_gradient(vars, grads);
  copy_to_gsl_vector(df, grads);
}

double GSLAdapter::my_f(const gsl_vector* v, void* params)
{
  auto this_adapter = reinterpret_cast<GSLAdapter*>(params);
  return this_adapter->my_f(v);
}

void GSLAdapter::my_df(const gsl_vector* v, void* params, gsl_vector* df)
{
  auto this_adapter = reinterpret_cast<GSLAdapter*>(params);
  this_adapter->my_df(v, df);
}

void GSLAdapter::my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df)
{
  auto this_adapter = reinterpret_cast<GSLAdapter*>(params);
  this_adapter->my_fdf(f, x, df);
}

void GSLAdapter::copy_to_gsl_vector(gsl_vector* gslv, const Eigen::VectorXd& eigenv)
{
  if (gslv->size != eigenv.size())
    throw std::runtime_error(
       fmt::format("copy_to_gsl_vector: vector size simatch {}!={}", gslv->size, eigenv.size()));
  for (size_t i = 0; i < eigenv.size(); ++i)
    gsl_vector_set(gslv, i, eigenv[i]);
}

void GSLAdapter::copy_from_gsl_vector(Eigen::VectorXd& eigenv, const gsl_vector* gslv)
{
  if (gslv->size != eigenv.size())
    throw std::runtime_error(
        fmt::format("copy_from_gsl_vector: vector size simatch {}!={}", eigenv.size(), gslv->size));
  for (size_t i = 0; i < eigenv.size(); ++i)
    eigenv[i] = gsl_vector_get(gslv, i);
}

}

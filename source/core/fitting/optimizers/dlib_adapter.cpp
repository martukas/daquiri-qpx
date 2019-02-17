#include <core/fitting/optimizers/dlib_adapter.h>
#include <dlib/matrix/matrix_mat.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

FitResult DLib::minimize(FittableFunction* fittable, double tolf)
{
  function_ = fittable;

  fitter_vector starting_point = dlib::mat(function_->variables());
  const auto& fe = std::bind(&DLib::eval, this, std::placeholders::_1);
  const auto& fd = std::bind(&DLib::derivative, this, std::placeholders::_1);

  dlib::find_min(dlib::bfgs_search_strategy(),
                 dlib::objective_delta_stop_strategy(1e-7).be_verbose(),
                 fe, fd, starting_point, -1);

  FitResult ret;

  ret.variables.setConstant(starting_point.size(), 0.0);
  for (long i = 0; i < starting_point.size(); ++i)
    ret.variables[i] = starting_point(i);

  return ret;
}

double DLib::eval(const DLib::fitter_vector& m) const
{
  Eigen::VectorXd v;
  v.setConstant(m.size(), 0.0);
  for (long i = 0; i < m.size(); ++i)
    v[i] = m(i);
  return function_->chi_sq(v);
}

DLib::fitter_vector DLib::derivative(const DLib::fitter_vector& m) const
{
  Eigen::VectorXd v;
  v.setConstant(m.size(), 0.0);
  for (long i = 0; i < m.size(); ++i)
    v[i] = m(i);
  Eigen::VectorXd g;
  g.setConstant(v.size(), 0.0);
  function_->chi_sq_gradient(v, g);
  return dlib::mat(g);
}

DLib::fitter_matrix DLib::hessian(const DLib::fitter_vector& m) const
{
  // \todo implement this
  (void) m;
  return fitter_matrix();
}

}

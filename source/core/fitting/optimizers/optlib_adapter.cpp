#include <core/fitting/optimizers/optlib_adapter.h>

#include <core/fitting/optimizers/cppoptlib/meta.h>
#include <core/fitting/optimizers/cppoptlib/problem.h>
#include <core/fitting/optimizers/cppoptlib/solver/bfgssolver.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

class OptlibFittableWrapper : public cppoptlib::Problem<double>
{
 public:
  FittableFunction* function_; /// < pointer to fittable function for function binding
  std::atomic<bool>* cancel_;

  using typename cppoptlib::Problem<double>::Scalar;
  using typename cppoptlib::Problem<double>::TVector;

  double value(const TVector& x) override
  {
    return function_->chi_sq(x);
  }

  void gradient(const TVector& x, TVector& grad) override
  {
    function_->chi_sq_gradient(x, grad);
  }

  bool callback(const cppoptlib::Criteria<double>& state, const TVector& x)
  {
    (void) state; // unused
    (void) x;     // unused
    return !cancel_->load();
  }
};

bool OptlibOptimizer::check_gradient(FittableFunction* fittable) const
{
  Eigen::VectorXd x = fittable->variables();

  OptlibFittableWrapper f;
  f.function_ = fittable;
  return f.checkGradient(x);
}

FitResult OptlibOptimizer::minimize(FittableFunction* fittable)
{
  OptlibFittableWrapper f;
  f.function_ = fittable;
  f.cancel_ = &cancel;
  Eigen::VectorXd x = fittable->variables();

  cppoptlib::BfgsSolver<OptlibFittableWrapper> solver;
  if (verbose)
    solver.setDebug(cppoptlib::DebugLevel::Low);
  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
  crit.iterations = maximum_iterations;
  crit.gradNorm = tolerance;
  solver.setStopCriteria(crit);

  solver.minimize(f, x);
  auto status = solver.status();

  FitResult ret;
  ret.converged = (status == cppoptlib::Status::GradNormTolerance)
      || (status == cppoptlib::Status::FDeltaTolerance)
      || (status == cppoptlib::Status::XDeltaTolerance);
  if (!ret.converged)
  {
    std::stringstream ss;
    ss << status;
    ret.error_message = ss.str();
  }
  ret.iterations = solver.criteria().iterations;
  f.finiteHessian(x, ret.inv_hessian);
  // \todo do we need to invert? normalize?
//  ret.inv_hessian = ret.inv_hessian.inverse();
  ret.variables = x;
  ret.value = fittable->chi_sq(x);
  return ret;
}

}

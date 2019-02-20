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

  bool callback(const cppoptlib::Criteria<double> &state, const TVector &x)
  {
    return !cancel_->load();
  }
};

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

  if (verbose)
  {
    std::stringstream ss;
    ss << solver.status();
    INFO("Optimization stopped with {}", ss.str());
  }

  FitResult ret;
  solver.minimize(f, x);
  ret.variables = x;
  ret.value = fittable->chi_sq(x);
  ret.converged = (solver.status() == cppoptlib::Status::GradNormTolerance)
      || (solver.status() == cppoptlib::Status::FDeltaTolerance)
      || (solver.status() == cppoptlib::Status::XDeltaTolerance);
  ret.iterations = solver.criteria().iterations;
  return ret;
}

}

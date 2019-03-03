#include <core/fitting/optimizers/optlib_adapter.h>

#include <core/fitting/optimizers/cppoptlib/meta.h>
#include <core/fitting/optimizers/cppoptlib/problem.h>

#include <core/fitting/optimizers/cppoptlib/solver/bfgssolver.h>
#include <core/fitting/optimizers/cppoptlib/solver/cmaessolver.h>
#include <core/fitting/optimizers/cppoptlib/solver/neldermeadsolver.h>

#include <core/util/custom_logger.h>

namespace DAQuiri
{

class OptlibFittableWrapper : public cppoptlib::Problem<double>
{
 public:
  OptlibFittableWrapper() = default;

  FittableFunction* function_; /// < pointer to fittable function for function binding
  std::atomic<bool>* cancel_;
  bool use_finite_gradient_{false};

  double epsilon{1e-10};
  double tolerance{1e-7};


  using typename cppoptlib::Problem<double>::Scalar;
  using typename cppoptlib::Problem<double>::TVector;

  double value(const TVector& x) override
  {
    return function_->chi_sq(x);
  }

  void gradient(const TVector& x, TVector& grad) override
  {
    if (!use_finite_gradient_)
      function_->chi_sq_gradient(x, grad);
    else
      finiteGradient(x, grad);
  }

  bool callback(const cppoptlib::Criteria<double>& state, const TVector& x) override
  {
    bool converged = false;
    if (std::isfinite(state.fDelta))
    {
      auto f = value(x);
      auto old_f = value(x) - state.fDelta;
      converged = ((2.0 * std::abs(state.fDelta)) <=
          (tolerance * (std::abs(f) + std::abs(old_f) + epsilon)));
//      INFO("{} = {} <= {}", converged, (2.0 * std::abs(state.fDelta)),
//           (tolerance * (std::abs(f) + std::abs(old_f) + epsilon)));
    }
    return !converged && !cancel_->load();
  }
};

bool OptlibOptimizer::check_gradient(FittableFunction* fittable) const
{
  Eigen::VectorXd x = fittable->variables();

  OptlibFittableWrapper f;
  f.function_ = fittable;

  return f.checkGradient(x);
}

void OptlibOptimizer::finite_gradient(FittableFunction* fittable,
                                      const Eigen::VectorXd& x,
                                      Eigen::VectorXd& gradients) const
{
  OptlibFittableWrapper f;
  f.function_ = fittable;
  f.use_finite_gradient_ = true;

  gradients.setConstant(x.size(), 0.0);
  f.gradient(x, gradients);
}


FitResult extract_status(const cppoptlib::Status& status, std::atomic<bool>* cancel)
{
  bool interrupted = cancel->load();
  FitResult ret;
  ret.converged = (status == cppoptlib::Status::GradNormTolerance)
      || (status == cppoptlib::Status::FDeltaTolerance)
      || (status == cppoptlib::Status::XDeltaTolerance)
      || ((status == cppoptlib::Status::UserDefined) && !interrupted);

  if (interrupted)
    ret.error_message = "Externally interrupted";
  else
  {
    std::stringstream ss;
    ss << status;
    ret.error_message = ss.str();
  }

  return ret;
}

FitResult solve(cppoptlib::BfgsSolver<OptlibFittableWrapper> solver,
                OptlibFittableWrapper f,
                OptlibOptimizer::GradientSelection grad_select);

FitResult solve(cppoptlib::BfgsSolver<OptlibFittableWrapper> solver,
                OptlibFittableWrapper f,
                OptlibOptimizer::GradientSelection grad_select)
{
  f.use_finite_gradient_ =
      (grad_select == OptlibOptimizer::GradientSelection::FiniteAlways);

  Eigen::VectorXd x = f.function_->variables();
  solver.minimize(f, x);

  FitResult ret = extract_status(solver.status(), f.cancel_);

  ret.iterations = solver.criteria().iterations;
  f.finiteHessian(x, ret.inv_hessian);
// \todo do we need to invert? normalize?
//  ret.inv_hessian = ret.inv_hessian.inverse();
  ret.variables = x;
  ret.value = f.function_->chi_sq(x);

  if (!ret.converged && (grad_select == OptlibOptimizer::GradientSelection::DefaultToFinite))
  {
    f.use_finite_gradient_ = true;
    ret = solve(solver, f, OptlibOptimizer::GradientSelection::FiniteAlways);
    ret.error_message = "Retry with finite gradient: " + ret.error_message;
    ret.used_finite_grads = true;
    f.use_finite_gradient_ = false;
  }

  return ret;
}

FitResult OptlibOptimizer::minimize(FittableFunction* fittable)
{
  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());

  cppoptlib::BfgsSolver<OptlibFittableWrapper> solver;
  if (verbose)
    solver.setDebug(cppoptlib::DebugLevel::Low);
  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
  crit.iterations = maximum_iterations;
  //crit.gradNorm = tolerance;
  crit.gradNorm = 0;
// \todo use different tolerance for this?
//  crit.xDelta = 1e-12;
  solver.setStopCriteria(crit);

  OptlibFittableWrapper f;
  f.function_ = fittable;
  f.cancel_ = &cancel;
  f.tolerance = tolerance;
  f.epsilon = epsilon;

  auto ret = solve(solver, f, gradient_selection);
  fittable->save_fit(ret);

  size_t perturbations = 0;
  while ((!ret.converged ||
          (perform_sanity_checks && !fittable->sane())) &&
      (perturbations < maximum_perturbations) &&
      fittable->perturb(random_generator))
  {
    perturbations++;
    ret = solve(solver, f, gradient_selection);
    fittable->save_fit(ret);
  }

  ret.perturbations = perturbations;

  return ret;
}

}

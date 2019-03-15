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
  OptlibFittableWrapper() = default;

  FittableFunction* function_; /// < pointer to fittable function for function binding
  std::atomic<bool>* cancel_;
  bool use_finite_gradient_{false};

  bool check_condition_ {false};
  double condition_tolerance_;
  double condition_epsilon_;

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
    if (check_condition_ && std::isfinite(state.fDelta))
    {
      auto f = value(x);
      auto old_f = value(x) - state.fDelta;
      converged = ((2.0 * std::abs(state.fDelta)) <=
          (condition_tolerance_ * (std::abs(f) + std::abs(old_f) + condition_epsilon_)));
//      INFO("{} = {} <= {}", converged, (2.0 * std::abs(state.fDelta)),
//           (tolerance * (std::abs(f) + std::abs(old_f) + epsilon)));
    }
    return !converged && !cancel_->load();
  }
};

std::string OptlibOptimizer::print_config(std::string prepend) const
{
  std::string ret;
  ret += prepend + fmt::format("verbosity level: {}\n", verbosity);
  ret += prepend + fmt::format("max iterations: {}\n", maximum_iterations);
  ret += prepend + fmt::format("min \u0394x: {}\n", min_x_delta);
  ret += prepend + fmt::format("min \u0394f: {}\n", min_f_delta);
  ret += prepend + fmt::format("min ||g||: {}\n", min_g_norm);
  ret += prepend + fmt::format("max Hessian condition: {}\n", max_condition);
  ret += prepend + fmt::format("use epsilon check: {}\n", (use_epsilon_check ? "YES" : "no"));
  if (use_epsilon_check)
  {
    ret += prepend + fmt::format("  tolerance: {}\n", tolerance);
    ret += prepend + fmt::format("  epsilon: {}\n", epsilon);
  }
  ret += prepend + fmt::format("max perturbations: {}\n", maximum_perturbations);
  ret += prepend + fmt::format("perform sanity checks: {}\n",
                               (perform_sanity_checks ? "YES" : "no"));
  return ret;
}

bool OptlibOptimizer::check_gradient(FittableFunction* fittable) const
{
  return check_gradient(fittable, fittable->variables());
}

bool OptlibOptimizer::check_gradient(FittableFunction* fittable,
                                     const Eigen::VectorXd& x) const
{
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
//      || (status == cppoptlib::Status::Condition)
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
                OptlibFittableWrapper& f, bool verbose,
                OptlibOptimizer::GradientSelection grad_select);

FitResult solve(cppoptlib::BfgsSolver<OptlibFittableWrapper> solver,
                OptlibFittableWrapper& f, bool verbose,
                OptlibOptimizer::GradientSelection grad_select)
{
  f.use_finite_gradient_ =
      (grad_select == OptlibOptimizer::GradientSelection::FiniteAlways);

  Eigen::VectorXd x = f.function_->variables();
  std::stringstream ss;
  solver.minimize(f, x, &ss);

  FitResult ret = extract_status(solver.status(), f.cancel_);
  ret.log = ss.str();

  ret.total_iterations = ret.iterations = solver.criteria().iterations;
  f.finiteHessian(x, ret.inv_hessian);
// \todo do we need to invert? normalize?
//  ret.inv_hessian = ret.inv_hessian.inverse();
  ret.variables = x;
  ret.value = f.function_->chi_sq(x);
  if (grad_select == OptlibOptimizer::GradientSelection::FiniteAlways)
    ret.total_finite_attempts++;
  else
    ret.total_analytic_attempts++;
  if (!ret.converged)
    ret.total_nonconvergences++;

  if (!ret.converged &&
      (grad_select == OptlibOptimizer::GradientSelection::DefaultToFinite))
  {
    auto ret_old = ret;
    if (verbose)
      ret_old.log +=
          "Retry with finite grad after failure to converge: " + ret.to_string(false) + "\n";
    ret = solve(solver, f, verbose, OptlibOptimizer::GradientSelection::FiniteAlways);
    ret.error_message = "Retry with finite gradient: " + ret.error_message;
    ret.used_finite_grads = true;
    ret.total_iterations += ret_old.total_iterations;
    ret.total_analytic_attempts += ret_old.total_analytic_attempts;
    ret.total_finite_attempts += ret_old.total_finite_attempts;
    ret.total_nonconvergences += ret_old.total_nonconvergences;
    ret.log = ret_old.log + ret.log;
  }

  return ret;
}

FitResult OptlibOptimizer::minimize(FittableFunction* fittable)
{
  std::mt19937 random_generator;
  random_generator.seed(std::random_device()());

  cppoptlib::BfgsSolver<OptlibFittableWrapper> solver;
  if (verbosity >= 5)
    solver.setDebug(cppoptlib::DebugLevel::High);
  else if (verbosity >= 4)
    solver.setDebug(cppoptlib::DebugLevel::Low);
  cppoptlib::Criteria<double> crit = cppoptlib::Criteria<double>::defaults();
  crit.iterations = maximum_iterations;
  crit.xDelta = min_x_delta;
  crit.fDelta = min_f_delta;
  crit.gradNorm = min_g_norm;
  crit.condition = max_condition;
  solver.setStopCriteria(crit);

  OptlibFittableWrapper f;
  f.function_ = fittable;
  f.cancel_ = &cancel;
  f.check_condition_ = use_epsilon_check;
  f.condition_tolerance_ = tolerance;
  f.condition_epsilon_ = epsilon;

  bool retry{false};
  size_t perturbations = 0;
  FitResult ret;
  do
  {
    auto old_ret = ret;
    retry = false;
    ret = solve(solver, f, (verbosity >= 2), gradient_selection);
    fittable->save_fit(ret);

    bool failed{false};
    if (!ret.converged)
    {
      // \todo print function itself
      if (verbosity >= 2)
        old_ret.log += "Failed to converge: " + ret.to_string(verbosity >= 6) + "\n";
      failed = true;
    }
    else if (perform_sanity_checks && !fittable->sane())
    {
      if (verbosity >= 2)
        old_ret.log += "Sanity check failed on: " + ret.to_string(verbosity >= 6) + "\n";
      ret.total_insane++;
      failed = true;
    }

    if (failed && (perturbations < maximum_perturbations))
    {
      retry = fittable->perturb(random_generator);
      if (verbosity >= 1)
      {
        // \todo print function itself
        if (!retry)
          old_ret.log += "Perturbation failed\n";
        else
        {
          std::stringstream ss;
          ss << f.function_->variables().transpose();
          old_ret.log += "Perturbed as:" + ss.str() + "\n";
        }
      }
      if (retry)
        perturbations++;
    }

    ret.total_iterations += old_ret.total_iterations;
    ret.total_analytic_attempts += old_ret.total_analytic_attempts;
    ret.total_finite_attempts += old_ret.total_finite_attempts;
    ret.total_nonconvergences += old_ret.total_nonconvergences;
    ret.total_insane += old_ret.total_insane;
    ret.log = old_ret.log + ret.log;
  }
  while (retry);

  ret.total_perturbations = perturbations;
  if (verbosity >= 1)
    ret.log = print_config() + ret.log;

  return ret;
}

}

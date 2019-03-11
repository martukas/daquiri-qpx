#pragma once

#include <atomic>
#include <core/fitting/optimizers/fittable_function.h>

namespace DAQuiri
{

/// \class AbstractOptimizer abstract_optimizer.h <core/fitting/optimizers/abstract_optimizer.h>
/// \brief Interface for an optimizer that can minimize an objective function
class AbstractOptimizer
{
 public:
  enum class GradientSelection
  {
    AnalyticalAlways,
    FiniteAlways,
    DefaultToFinite
  } gradient_selection;

  std::atomic<bool> cancel{false}; /// < to stop the optimizer from any thread

  size_t verbosity{0};

  size_t maximum_iterations{3000};

  double min_x_delta{0.};
  double min_f_delta{0.};
  double min_g_norm{0.};
  double max_condition{0.};

  bool use_epsilon_check{true};
  double tolerance{1e-10};
  double epsilon{100 * std::numeric_limits<double>::min()};

  size_t maximum_perturbations{10};
  bool perform_sanity_checks{false};

  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable instance of an objective FittableFunction to be minimized
  virtual FitResult minimize(FittableFunction* fittable) = 0;

  /// \brief calculates finite gradient for supplied function
  /// \param fittable instance of an objective FittableFunction
  /// \param x function variables at which to evaluate
  /// \param gradients output vector for gradients
  virtual void finite_gradient(FittableFunction* fittable,
                               const Eigen::VectorXd& x,
                               Eigen::VectorXd& gradients) const
  {
    (void) fittable; // < unused
    (void) x; // < unused
    (void) gradients; // < unused
  }

  /// \brief checks if analytical gradient is ok
  /// \returns true if gradiengt is probably ok, not a guarantee
  /// \param fittable instance of an objective FittableFunction
  virtual bool check_gradient(FittableFunction* fittable) const
  {
    (void) fittable; // < unused
    return false;
  }

  /// \brief checks if analytical gradient is ok
  /// \returns true if gradiengt is probably ok, not a guarantee
  /// \param fittable instance of an objective FittableFunction
  /// \param x function variables at which to evaluate
  virtual bool check_gradient(FittableFunction* fittable, const Eigen::VectorXd& x) const
  {
    (void) fittable; // < unused
    (void) x; // < unused
    return false;
  }
};

}

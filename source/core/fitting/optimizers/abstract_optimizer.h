#pragma once

#include <atomic>
#include <core/fitting/optimizers/fittable_function.h>
#include <core/fitting/optimizers/fit_result.h>

namespace DAQuiri
{

/// \class AbstractOptimizer abstract_optimizer.h <core/fitting/optimizers/abstract_optimizer.h>
/// \brief Interface for an optimizer that can minimize an objective function
class AbstractOptimizer
{
 public:
  std::atomic<bool> cancel{false}; /// < to stop the optimizer from any thread

  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  /// \param tolf ???
  virtual FitResult minimize(FittableFunction* fittable, double tolf) = 0;
};

}

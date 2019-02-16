#pragma once

#include <atomic>
#include <core/fitting/optimizers/fittable_function.h>
#include <core/fitting/optimizers/fit_results.h>

namespace DAQuiri
{

class AbstractOptimizer
{
 public:
  std::atomic<bool> cancel{false};
  virtual FitResult minimize(FittableFunction* fittable, double tolf) = 0;
};

}

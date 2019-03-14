#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri
{

// \todo add locks

/// \class BFGS BFGS.h <core/fitting/optimizers/BFGS.h>
/// \brief Implementation of the full-memory Broyden-Fletcher-Goldfarb-Shanno optimizer.
class BudapestOptimizer : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  FitResult minimize(FittableFunction* fittable) override;

  double eps{1e-10};
};

}

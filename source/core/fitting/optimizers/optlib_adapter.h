#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri
{

// \todo add locks

/// \class OptlibOptimizer dlib_adapter.h <core/fitting/optimizers/dlib_adapter.h>
/// \brief Interface adapter for the cppoptlib implementation of BFGS.
///         Does not derive uncertainties.
class OptlibOptimizer : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  FitResult minimize(FittableFunction* fittable) override;

  bool check_gradient(FittableFunction* fittable) const;

  enum class GradientSelection
  {
    AnalyticalAlways,
    FiniteAlways,
    DefaultToFinite
  };

  GradientSelection gradient_selection;
  bool perform_sanity_checks {false};

 private:


};

}

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
  /// \param fittable instance of an objective FittableFunction to be minimized
  FitResult minimize(FittableFunction* fittable) override;

  /// \brief calculates finite gradient for supplied function
  /// \param fittable instance of an objective FittableFunction
  /// \param x function variables at which to evaluate
  /// \param gradients output vector for gradients
  void finite_gradient(FittableFunction* fittable,
                       const Eigen::VectorXd& x,
                       Eigen::VectorXd& gradients) const override;

  /// \brief checks if analytical gradient is ok
  /// \returns true if gradiengt is probably ok, not a guarantee
  /// \param fittable instance of an objective FittableFunction
  bool check_gradient(FittableFunction* fittable) const override;

  /// \brief checks if analytical gradient is ok
  /// \returns true if gradiengt is probably ok, not a guarantee
  /// \param fittable instance of an objective FittableFunction
  bool check_gradient(FittableFunction* fittable,
                      const Eigen::VectorXd& x) const override;

  std::string print_config(std::string prepend = "") const;
};

}

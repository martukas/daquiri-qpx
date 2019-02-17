#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

#include <dlib/optimization.h>

namespace DAQuiri
{

// \todo add locks

/// \class DLibOptimizer dlib_adapter.h <core/fitting/optimizers/dlib_adapter.h>
/// \brief Interface adapter for the dlib implementation of BFGS.
///         Does not derive uncertainties.
class DLibOptimizer : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  FitResult minimize(FittableFunction* fittable) override;

  double minimum_value{-1.0};

 private:
  FittableFunction* function_; /// < pointer to fittable function for function binding
  Eigen::VectorXd variables_;
  Eigen::VectorXd gradients_;

  using fitter_vector = dlib::matrix<double, 0, 1>;

  double eval(const fitter_vector& vars);
  fitter_vector derivative(const fitter_vector& vars);

};

}

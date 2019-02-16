#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

#include <dlib/optimization.h>


namespace DAQuiri
{

// \todo add locks

/// \class DLib dlib_adapter.h <core/fitting/optimizers/dlib_adapter.h>
/// \brief Interface adapter for the dlib implementation of BFGS. Does not derive uncertainties.
class DLib : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  /// \param tolf ???
  FitResult minimize(FittableFunction* fittable, double tolf) override;

 private:
  FittableFunction* function_; /// < pointer to fittable function for function binding

  using fitter_vector = dlib::matrix<double, 0, 1>;
  using fitter_matrix = dlib::matrix<double>;

  double eval(const fitter_vector& m) const;
  fitter_vector derivative (const fitter_vector& m) const;
  fitter_matrix hessian(const fitter_vector& m) const;

};

}

#pragma once

#pragma GCC diagnostic push
#ifdef __GNUC__
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

#include <random>

namespace DAQuiri
{

/// \class FittableFunction fittable_function.h <core/fitting/optimizers/fittable_function.h>
/// \brief Interface for an objective function that can be minimized by the optimizer.
class FittableFunction
{
 public:
  FittableFunction() = default;
  virtual ~FittableFunction() = default;

  /// \returns a vectorized list of all variables in the model. Only variables considered
  ///          for optimization should be included.
  virtual Eigen::VectorXd variables() const = 0;

  /// \brief Evaluates the objective function with the provided set of variables.
  /// \returns chi-squared of the fit as compared to the empirical data.
  /// \param fit a vector containing variable values being tested.
  virtual double chi_sq(const Eigen::VectorXd& fit) const = 0;

  /// \brief Evaluates the objective function with the provided set of variables and
  ///         simultaneously provides gradients for each variable at this point.
  /// \returns chi-squared of the fit as compared to the empirical data.
  /// \param fit a vector containing variable values being tested.
  /// \param gradients a vector into which gradients for each variable will be written.
  virtual double chi_sq_gradient(const Eigen::VectorXd& fit,
                                 Eigen::VectorXd& gradients) const = 0;

  // \todo document this
  virtual bool perturb(std::mt19937&) { return false; }

//  virtual double chi_sq_gradient_hessian(const column_vector& x, column_vector& der,
//      general_matrix& hess) const = 0;
};

}

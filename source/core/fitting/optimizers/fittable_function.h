#pragma once

#include <core/util/eigen_fix.h>
#include <random>
#include <core/fitting/optimizers/fit_result.h>


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
  virtual double chi_sq(const Eigen::VectorXd& fit) const
  {
    // \todo add notes for implementers interested in performance
    // \todo only resize if necessary
    dummy_gradient.setConstant(fit.size(), 0.);
    return this->chi_sq_gradient(fit, dummy_gradient);
  }

  /// \brief Evaluates the objective function with the provided set of variables and
  ///         simultaneously provides gradients for each variable at this point.
  /// \returns chi-squared of the fit as compared to the empirical data.
  /// \param fit a vector containing variable values being tested.
  /// \param gradients a vector into which gradients for each variable will be written.
  virtual double chi_sq_gradient(const Eigen::VectorXd& fit,
                                 Eigen::VectorXd& gradients) const = 0;

  // \todo make this function pure virtual?
  /// \brief saves fit result, possibly with uncertainties
  /// \param result result of optimization attempt
  virtual void save_fit(const FitResult&) {};

  // \todo document this
  virtual bool perturb(std::mt19937&) { return false; }

  // \todo document this
  virtual bool sane() const { return true; }

//  virtual double chi_sq_gradient_hessian(const column_vector& x, column_vector& der,
//      general_matrix& hess) const = 0;
 protected:
  mutable Eigen::VectorXd dummy_gradient;
};

}

#pragma once

#include <core/fitting/optimizers/fittable_function.h>
#include <core/fitting/optimizers/fit_result.h>
#include <core/fitting/weighted_data.h>

namespace DAQuiri
{

/// \class FittableRegion fittable_region.h <core/fitting/fittable_region.h>
/// \brief Partially implements objective FittableFunction for a spectrum region, and provides
///         interface for supplying a model function to be compared against the data.
///         Five new virtual functions must be implemented as well as the variables() function
///         required by FittableFunction.
class FittableRegion : public FittableFunction
{
 public:
  int32_t variable_count {0}; /// < counter for indexing/mapping of model variables to vectors
  WeightedData data;          /// < empirical data

  FittableRegion() = default;
  virtual ~FittableRegion() = default;

  /// \brief updates variable indices for vector mapping, updates variable_count
  virtual void update_indices() = 0;

  /// \brief evaluates the current state of the model function at a single channel
  /// \returns value of function at channel
  /// \param chan channel value
  virtual double eval(double chan) const = 0;

  /// \brief evaluates a particular fit of the model function at a single channel
  /// \returns value of function at channel
  /// \param chan channel value
  /// \param fit vector of model function variables to be evaluated
  virtual double eval_at(double chan, const Eigen::VectorXd& fit) const
  {
    // \todo add notes for implementers interested in performance
    // \todo only resize if necessary
    dummy_gradient.setConstant(fit.size(), 0.);
    return this->eval_grad_at(chan, fit, dummy_gradient);
  }

  /// \brief evaluates a particular fit of the model at a single channel and calculates
  ///        the function gradients at the same channel
  /// \returns value of function at channel
  /// \param chan channel value
  /// \param fit vector of model function variables to be evaluated
  /// \param grads vector into which channel gradients will be written
  virtual double eval_grad_at(double chan, const Eigen::VectorXd& fit,
      Eigen::VectorXd& grads) const = 0;

  // \todo document this
  virtual std::string to_string(std::string prepend = "") const = 0;

  /// \returns chi squared of the current state of the model function as compared to
  ///             the empirical data
  double chi_sq() const;

  /// \returns chi squared of of a particular fit of the model function as compared to
  ///         the empirical data
  /// \param fit vector of model function variables to be evaluated
  double chi_sq(const Eigen::VectorXd& fit) const override;

  /// \returns chi squared of of a particular fit of the model function as compared to
  ///         the empirical data
  /// \param fit vector of model function variables to be evaluated
  /// \param gradients vector into which variable gradients will be written
  double chi_sq_gradient(const Eigen::VectorXd& fit,
                         Eigen::VectorXd& gradients) const override;

  /// \returns degrees of freedom for the objective function being evaluated
  /// \pre number of variables in the model function must be less than number of channels in region
  double degrees_of_freedom() const;
};

}
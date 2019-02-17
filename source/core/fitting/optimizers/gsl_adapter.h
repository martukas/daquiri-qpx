#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

#include <gsl/gsl_multimin.h>

namespace DAQuiri
{

// \todo add locks

/// \class DLib dlib_adapter.h <core/fitting/optimizers/dlib_adapter.h>
/// \brief Interface adapter for the dlib implementation of BFGS. Does not derive uncertainties.
class GSLAdapter : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  /// \param tolf ???
  FitResult minimize(FittableFunction* fittable, double tolf) override;

 private:
  FittableFunction* function_; /// < pointer to fittable function for function binding

  Eigen::VectorXd vars;
  Eigen::VectorXd grads;

  double my_f(const gsl_vector* v);
  void my_df(const gsl_vector* v, gsl_vector* df);
  void my_fdf(double* ret, const gsl_vector* x, gsl_vector* df);

  static double my_f(const gsl_vector* v, void* params);
  static void my_df(const gsl_vector* v, void* params, gsl_vector* df);
  static void my_fdf(const gsl_vector* x, void* params, double* f, gsl_vector* df);

  static void copy_to_gsl_vector(gsl_vector* gslv, const Eigen::VectorXd& eigenv);
  static void copy_from_gsl_vector(Eigen::VectorXd& eigenv, const gsl_vector* gslv);
};

}

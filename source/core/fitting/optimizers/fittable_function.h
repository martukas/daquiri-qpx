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


namespace DAQuiri
{

class FittableFunction
{
 public:
  FittableFunction() = default;
  virtual ~FittableFunction() = default;

  virtual Eigen::VectorXd variables() const = 0;
  virtual double chi_sq(const Eigen::VectorXd& fit) const = 0;
  virtual double chi_sq_gradient(const Eigen::VectorXd& fit,
                                 Eigen::VectorXd& gradients) const = 0;

//  virtual double chi_sq_gradient_hessian(const column_vector& x, column_vector& der,
//      general_matrix& hess) const = 0;
};

}

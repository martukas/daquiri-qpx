#pragma once

#include <cstdint>
#include <vector>

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

struct FitResult
{
  Eigen::VectorXd variables;
  Eigen::SparseMatrix<double> inv_hessian;
  size_t iterations{0};
  bool converged{false};
};

class Fittable
{
 public:
  Fittable() = default;
  virtual ~Fittable() = default;

  virtual double operator()(const Eigen::VectorXd& fit,
                            Eigen::VectorXd& gradients) const = 0;

  virtual Eigen::VectorXd variables() const = 0;
  virtual double chi_sq(const Eigen::VectorXd& fit) const = 0;

//  void get_derivative_and_hessian (const column_vector& x, column_vector& der,
//      general_matrix& hess) const
};

}

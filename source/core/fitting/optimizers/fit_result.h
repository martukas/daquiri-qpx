#pragma once

#include <core/util/eigen_fix.h>
#include <string>

namespace DAQuiri
{

using hessian_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

/// \struct FitResult fit_result.h <core/fitting/optimizers/fit_result.h>
/// \brief provides all relevant information from the results of an optimization attempt.
struct FitResult
{
  // \todo check to confirm all variables are finite

  Eigen::VectorXd variables;       /// < variable values arrived at
  hessian_t inv_hessian;           /// < inverse Hessian matrix of recent fit iterations
  bool converged{false};           /// < whether convergence was achieved
  size_t iterations{0};            /// < number of iterations used to reach result
  double value;                    /// < most recent evaluation result
  bool used_finite_grads {false};

  std::string error_message;
  size_t total_perturbations{0};
  size_t total_iterations{0};
  size_t total_analytic_attempts{0};
  size_t total_finite_attempts{0};
  size_t total_nonconvergences{0};
  size_t total_insane{0};
  //double time_elapsed{0.0};

  std::string log;

  std::string to_string(bool with_hessian = false) const;
};

}

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

#include <string>

namespace DAQuiri
{

/// \struct FitResult fit_result.h <core/fitting/optimizers/fit_result.h>
/// \brief provides all relevant information from the results of an optimization attempt.
struct FitResult
{
  Eigen::VectorXd variables;               /// < variable values arrived at
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
      inv_hessian; /// < inverse Hessian matrix of recent fit iterations
  bool converged{false};                   /// < whether convergence was achieved
  size_t iterations{0};                    /// < number of iterations used to reach result
  double value;                            /// < most recent evaluation result
  bool used_finite_grads {false};

  std::string error_message;
  size_t total_perturbations{0};
  size_t total_iterations{0};
  size_t total_analytic_attempts{0};
  size_t total_finite_attempts{0};
  size_t total_nonconvergences{0};
  size_t total_insane{0};
  //double time_elapsed{0.0};

  std::string to_string(bool with_hessian = false) const;
};

}

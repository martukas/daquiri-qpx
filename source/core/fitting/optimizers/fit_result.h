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
  size_t iterations{0};                    /// < number of iterations used to reach result
  bool converged{false};                   /// < whether convergence was achieved
  double value;                            /// < most recent evaluation result

  std::string to_string(bool with_hessian = false) const;
};

}

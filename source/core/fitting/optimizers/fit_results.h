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

struct FitResult
{
  Eigen::VectorXd variables;
  Eigen::SparseMatrix<double> inv_hessian;
  size_t iterations{0};
  bool converged{false};
};

}

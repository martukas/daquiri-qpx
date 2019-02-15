#pragma once

#include <atomic>
#include <core/fitting/fittable_function.h>

namespace DAQuiri
{

struct FitResult
{
  Eigen::VectorXd variables;
  Eigen::SparseMatrix<double> inv_hessian;
  size_t iterations{0};
  bool converged{false};
};

class AbstractOptimizer
{
 public:
  std::atomic<bool> cancel{false};
  virtual FitResult BFGSMin(Fittable* fittable, double tolf) = 0;
};

}

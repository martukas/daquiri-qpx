#pragma once

#include <cstdint>
#include <vector>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

namespace Hypermet
{

struct FitResult
{
  std::vector<double> variables;
  Eigen::SparseMatrix<double> inv_hessian;
  size_t iterations{0};
  bool converged{false};
};

class Fittable
{
 public:
  Fittable() = default;
  virtual ~Fittable() = default;

  virtual std::vector<double> variables() const = 0;
  virtual double degrees_of_freedom() const = 0;
  virtual double chi_sq(const std::vector<double>& fit) const = 0;
  virtual double grad_chi_sq(const std::vector<double>& fit,
                             std::vector<double>& gradients) const = 0;
};

}

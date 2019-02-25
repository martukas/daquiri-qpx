#include <core/fitting/optimizers/fit_result.h>
#include <core/util/UTF_extensions.h>
#include <fmt/format.h>

namespace DAQuiri
{

std::string FitResult::to_string(bool with_hessian) const
{
  std::stringstream ss;
  ss << variables.transpose();
  auto ret =
      fmt::format("{}{}{} at val={} after {} iterations  vars={} {}",
                  (converged ? "CONVERGED" : "NO CONVERGENCE"),
                  (used_finite_grads ? " FINITE" : ""),
                  ((perturbations > 0) ? (" Perturbation:" + std::to_string(perturbations)) : ""),
                  value, iterations, ss.str(), error_message);

  if (with_hessian && inv_hessian.innerSize() && inv_hessian.outerSize())
  {
    std::stringstream ss2;
    ss2 << "\n" << inv_hessian;
    ret += ss2.str();
  }

  return ret;
}

}

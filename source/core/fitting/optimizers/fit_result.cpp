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
      fmt::format("{}{}{} at f={} after {} iterations  x={} \"{}\"",
                  (converged ? "CONVERGED" : "NO CONVERGENCE"),
                  (used_finite_grads ? " FINITE" : ""),
                  ((total_perturbations > 0) ? (" Perturbation:" + std::to_string(total_perturbations)) : ""),
                  value, iterations, ss.str(), error_message);

  ret += fmt::format(" tot_iter={} analytical={} finite={} nonconv={} insane={}",
      total_iterations, total_analytic_attempts, total_finite_attempts, total_nonconvergences,
      total_insane);

  if (with_hessian && inv_hessian.innerSize() && inv_hessian.outerSize())
  {
    std::stringstream ss2;
    ss2 << "\n" << inv_hessian;
    ret += ss2.str();
  }

  return ret;
}

}

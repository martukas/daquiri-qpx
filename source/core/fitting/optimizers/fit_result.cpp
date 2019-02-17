#include <core/fitting/optimizers/fit_result.h>
#include <fmt/format.h>

namespace DAQuiri
{

std::string FitResult::to_string() const
{
  std::stringstream ss;
  ss << variables.transpose();
  return fmt::format("{} at {} iterations={}  vars={}",
                     (converged ? "CONVERGED" : "NO CONVERGENCE"),
                     value, iterations, ss.str());
}

}

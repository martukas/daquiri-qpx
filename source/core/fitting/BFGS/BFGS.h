#pragma once

#include <core/fitting/BFGS/Fittable.h>
#include <atomic>

namespace DAQuiri
{

// \todo add locks
class BFGS
{
 public:
  std::atomic<bool> cancel{false};
  FitResult BFGSMin(const Fittable* const fittable, double tolf);

 private:
  double Sign(double a, double b);
  double BrentDeriv(const Fittable* const fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const std::vector<double>& variables, const std::vector<double>& hessian);
  void Bracket(const Fittable* const fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const std::vector<double>& variables, const std::vector<double>& hessian);
  double fgv(const Fittable* const fittable, double lambda,
             std::vector<double> variables, std::vector<double> hessian);
  double dfgv(const Fittable* const fittable, double lambda,
              std::vector<double> variables, std::vector<double> hessian);
  double LinMin(const Fittable* const fittable,
                std::vector<double>& variables,
                std::vector<double> hessian);
};

}

#pragma once

#include <core/fitting/BFGS/Fittable.h>
#include <atomic>

namespace Hypermet
{

class BFGS
{
 public:
  FitResult BFGSMin(const Fittable *const fittable, double tolf, std::atomic<bool>& cancel);

 private:
  double Sign(double a, double b);
  double BrentDeriv(const Fittable *const fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const std::vector<double>& variables, const std::vector<double>& hessian);
  void Bracket(const Fittable *const fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const std::vector<double>& variables, const std::vector<double>& hessian);
  double fgv(const Fittable *const fittable, double lambda,
      std::vector<double> variables, std::vector<double> hessian);
  double dfgv(const Fittable *const fittable, double lambda,
      std::vector<double> variables, std::vector<double> hessian);
  double LinMin(const Fittable *const fittable,
                std::vector<double>& variables,
                std::vector<double> hessian);
};

}

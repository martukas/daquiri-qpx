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
  FitResult BFGSMin(Fittable* fittable, double tolf);

 private:
  double Sign(double a, double b);
  double BrentDeriv(Fittable* fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const Eigen::VectorXd& variables,
                    const Eigen::VectorXd& hessian,
                    Eigen::VectorXd& chan_gradients);
  void Bracket(Fittable* fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const Eigen::VectorXd& variables, const Eigen::VectorXd& hessian);
  double fgv(Fittable* fittable, double lambda,
             Eigen::VectorXd variables, Eigen::VectorXd hessian);
  double dfgv(Fittable* fittable, double lambda,
              Eigen::VectorXd variables, Eigen::VectorXd hessian,
              Eigen::VectorXd& chan_gradients);
  double LinMin(Fittable* fittable,
                Eigen::VectorXd& variables,
                Eigen::VectorXd hessian,
                Eigen::VectorXd& chan_gradients);
};

}

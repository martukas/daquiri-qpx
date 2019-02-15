#pragma once

#include <core/fitting/abstract_optimizer.h>

namespace DAQuiri
{

// \todo add locks
class BFGS : public AbstractOptimizer
{
 public:
  FitResult BFGSMin(Fittable* fittable, double tolf) override;

 private:
  double Sign(double a, double b);
  double BrentDeriv(Fittable* fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const Eigen::VectorXd& variables,
                    const Eigen::VectorXd& hessian);
  void Bracket(Fittable* fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const Eigen::VectorXd& variables, const Eigen::VectorXd& hessian);
  double fgv(Fittable* fittable, double lambda,
             Eigen::VectorXd variables, Eigen::VectorXd hessian);
  double dfgv(Fittable* fittable, double lambda,
              Eigen::VectorXd variables, Eigen::VectorXd hessian);
  double LinMin(Fittable* fittable,
                Eigen::VectorXd& variables,
                Eigen::VectorXd hessian);
};

}

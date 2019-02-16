#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri
{

// \todo add locks

/// \class BFGS BFGS.h <core/fitting/optimizers/BFGS.h>
/// \brief Implementation of the full-memory Broyden-Fletcher-Goldfarb-Shanno optimizer.
class BFGS : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  /// \param tolf ???
  FitResult minimize(FittableFunction* fittable, double tolf) override;

 private:
  double Sign(double a, double b);
  double BrentDeriv(FittableFunction* fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const Eigen::VectorXd& variables,
                    const Eigen::VectorXd& hessian);
  void Bracket(FittableFunction* fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const Eigen::VectorXd& variables, const Eigen::VectorXd& hessian);
  double fgv(FittableFunction* fittable, double lambda,
             Eigen::VectorXd variables, Eigen::VectorXd hessian);
  double dfgv(FittableFunction* fittable, double lambda,
              Eigen::VectorXd variables, Eigen::VectorXd hessian);
  double LinMin(FittableFunction* fittable,
                Eigen::VectorXd& variables,
                Eigen::VectorXd hessian);
};

}

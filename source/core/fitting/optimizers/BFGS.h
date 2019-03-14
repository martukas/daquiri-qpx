#pragma once

#include <core/fitting/optimizers/abstract_optimizer.h>

namespace DAQuiri
{

// \todo add locks

/// \class BFGS BFGS.h <core/fitting/optimizers/BFGS.h>
/// \brief Implementation of the full-memory Broyden-Fletcher-Goldfarb-Shanno optimizer.
class BudapestOptimizer : public AbstractOptimizer
{
 public:
  /// \brief minimizes a supplied function
  /// \returns result of the optimization attempt
  /// \param fittable a concrete instance of an objective FittableFunction to be minimized
  FitResult minimize(FittableFunction* fittable) override;

  double eps{1e-10};

  size_t brent_maximum_iterations{500};
  double brent_zeps{1e-10};

  double bracket_glimit{100.0};
  double bracket_tiny{1e-20};

  double linmin_tol{0.0001};

 private:
  double Sign(double a, double b);
  double BrentDeriv(FittableFunction* fittable,
                    double a, double b, double c, double tol, double& xmin,
                    const Eigen::VectorXd& variables,
                    const Eigen::VectorXd& search_direction);
  void Bracket(FittableFunction* fittable,
               double& a, double& b, double& c, double& fa, double& fb, double& fc,
               const Eigen::VectorXd& variables, const Eigen::VectorXd& search_direction);
  double fgv(FittableFunction* fittable, double lambda,
             Eigen::VectorXd variables, Eigen::VectorXd search_direction);
  double dfgv(FittableFunction* fittable, double lambda,
              Eigen::VectorXd variables, Eigen::VectorXd search_direction);
  double LinMin(FittableFunction* fittable,
                Eigen::VectorXd& variables,
                Eigen::VectorXd search_direction);
};

}

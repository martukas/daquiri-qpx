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

  double linmin_tolerance{0.0001};

 private:
  struct StepEval
  {
    StepEval(FittableFunction* fittable, const Eigen::VectorXd& variables, const Eigen::VectorXd& search_direction);

    StepEval(const StepEval& other) = default;

    StepEval& operator=(const StepEval& other);

    void recalc_f(double lambda);

    void recalc_df(double lambda);

    FittableFunction* fittable_;
    const Eigen::VectorXd* variables_;
    const Eigen::VectorXd* search_direction_;

    double step;
    double f;
    double dot;
  };

  double Sign(double a, double b);
  StepEval BrentDeriv(FittableFunction* fittable,
                      double lambda1,
                      double initial_lambda,
                      double lambda2,
                      const Eigen::VectorXd &variables,
                      const Eigen::VectorXd &search_direction);

  void Bracket(FittableFunction* fittable,
               StepEval& a_step, StepEval& b_step, StepEval& c_step,
               const Eigen::VectorXd& variables, const Eigen::VectorXd& search_direction);

  double LinMin(FittableFunction* fittable,
                Eigen::VectorXd& variables,
                const Eigen::VectorXd& search_direction);
};

}

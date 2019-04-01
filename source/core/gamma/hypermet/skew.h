#pragma once

#include <core/fitting/data_model/value.h>
#include <core/gamma/hypermet/gaussian_precalc_values.h>


namespace DAQuiri
{

class Tail
{
 public:
  Tail() = default;
  Tail(Side s) : side(s) {}

  bool override{false};
  bool enabled{true};
  Side side {Side::left};
  SineBoundedValue amplitude, slope;

  void reset_indices();

  /// \brief if enabled, saves and increments indices for variables flagged for fitting
  /// \param idx index from parent model
  void update_indices(int32_t& i);

  /// \brief writes variables into fit vector for optimization, only if indices are valid
  /// \param fit variables to be optimized
  void put(Eigen::VectorXd& fit) const;

  /// \brief reads variables from fit vector, only if indices are valid
  /// \param fit variables currently being evaluated by optimizer
  void get(const Eigen::VectorXd& fit);

  /// \brief calculates and stores uncertainties based on current fit
  /// \param diagonals of the inverse Hessian matrix of the current fit
  /// \param chisq_norm the normalized Chi-squared of the current fit
  void get_uncerts(const Eigen::VectorXd& diagonals, double chisq_norm);

  /// \returns current value of the tail function
  /// \param pre-calculated intermediate values based on related Gaussian function
  double eval(const PrecalcVals& pre) const;

  /// \returns value of the tail function at current fit
  /// \param pre-calculated intermediate values based on related Gaussian function
  /// \param fit variables currently being evaluated by optimizer
  double eval_at(const PrecalcVals& pre, const Eigen::VectorXd& fit) const;

  /// \returns current gradient (slope) of the tail function
  /// \param pre-calculated intermediate values based on related Gaussian function
  /// \param grads gradients for current fit
  double eval_grad(const PrecalcVals& pre,
                   Eigen::VectorXd& grads) const;

  /// \returns gradient (slope) of the tail function at current fit
  /// \param pre-calculated intermediate values based on related Gaussian function
  /// \param fit variables currently being evaluated by optimizer
  /// \param grads gradients for current fit
  double eval_grad_at(const PrecalcVals& pre, const Eigen::VectorXd& fit,
                      Eigen::VectorXd& grads) const;

  bool sane(double amp_min_epsilon, double amp_max_epsilon, double slope_epsilon) const;

  std::string to_string() const;

 private:
  double eval_with(const PrecalcVals& pre, double ampl, double slp) const;

};

void to_json(nlohmann::json& j, const Tail& s);
void from_json(const nlohmann::json& j, Tail& s);


}

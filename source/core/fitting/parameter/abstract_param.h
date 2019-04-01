#pragma once

#include <cinttypes>
#include <nlohmann/json.hpp>
#include <core/util/eigen_fix.h>

namespace DAQuiri
{

// \todo invalidate uncertainties automatically when changing value
// \todo sanity checks for nan and inf

class AbstractValue
{
 public:
  AbstractValue() = default;
  virtual ~AbstractValue() = default;

  bool to_fit{true}; ///< flags Value for fitting

  static constexpr int32_t InvalidIndex {-1};

  /// \brief if Value is flagged for fitting, saves and increments index
  ///        otherwise, invalidates internally stored index
  /// \param idx index from parent model
  void update_index(int32_t& idx);

  /// \brief invalidates internally stored index
  void reset_index();

  /// \returns internally stored index for fitting
  int32_t index() const;

  /// \returns true if index is valid
  bool valid_index() const;

  /// \returns current proxy variable
  double x() const;
  /// \brief sets current proxy variable
  void x(double new_x);

  /// \returns current nominal value
  double val() const;

  /// \returns value as function of current fit, if variable index is valid
  ///          else, returns current nominal value
  /// \param fit variables currently being evaluated by optimizer
  double val_from(const Eigen::VectorXd& fit) const;

  /// \brief sets current proxy variable so that nominal value equals new_val
  /// \param new_val new nominal value to set
  virtual void val(double new_val) = 0;

  /// \returns transforms proxy variable into nominal value
  /// \param at_x proxy variable to be evaluated
  virtual double val_at(double at_x) const = 0;

  /// \returns current gradient (slope), as function of currently stored proxy variable
  double grad() const;

  /// \returns gradient(slope) of variable in current fit, if variable index is valid
  ///          else, returns current gradient (slope)
  /// \param fit variables currently being evaluated by optimizer
  double grad_from(const Eigen::VectorXd& fit) const;

  /// \returns transforms proxy variable into gradient (slope)
  /// \param at_x proxy variable to be evaluated
  virtual double grad_at(double at_x) const = 0;

  /// \returns current uncertainty
  double uncert() const;

  /// \brief sets current uncertainty
  void uncert(double new_uncert);

  /// \brief writes current proxy variable into fit vector for optimization, only
  ///        if variable index is valid
  /// \param fit variables to be optimized
  void put(Eigen::VectorXd& fit) const;

  /// \brief reads current proxy variable from fit vector, only if variable index is valid
  /// \param fit variables currently being evaluated by optimizer
  void get(const Eigen::VectorXd& fit);

  /// \brief calculates and stores uncertainty for the value based on current fit
  /// \param diagonals of the inverse Hessian matrix of the current fit
  /// \param chisq_norm the normalized Chi-squared of the current fit
  void get_uncert(const Eigen::VectorXd& diagonals, double chisq_norm);

  virtual std::string to_string() const;
  friend void to_json(nlohmann::json& j, const AbstractValue& s);
  friend void from_json(const nlohmann::json& j, AbstractValue& s);

 private:
  double x_{0.0};

  // \todo what's this for?
  double dx_{0.0};

  double val_uncert_{0.0};
  int32_t index_{InvalidIndex};
};

}

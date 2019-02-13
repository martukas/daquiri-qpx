#pragma once

#include <cinttypes>
#include <vector>
#include <nlohmann/json.hpp>

#pragma GCC diagnostic push
#ifdef __GNUC__
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wmisleading-indentation"
#endif
#endif
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#include <Eigen/Sparse>
#pragma GCC diagnostic pop

namespace DAQuiri
{

class AbstractValue
{
 public:
  AbstractValue() = default;
  virtual ~AbstractValue() = default;

  bool to_fit{true}; ///< flags Value for fitting

  /// \brief if Value is flagged for fitting, saves and increments index
  ///        otherwise, invalidates internally stored index
  /// \param idx index from parent model
  void update_index(int32_t& idx);

  // \todo document this
  void reset_index();

  /// \returns internally stored index for fitting
  int32_t index() const;

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

  // \todo make this private
 protected:
  double x_{0.0};

  // \todo what's this for?
  double dx_{0.0};

  double val_uncert_{0.0};
  int32_t index_{-1};
};

class Value : public AbstractValue
{
 public:
  Value() = default;

  using AbstractValue::x;
  using AbstractValue::val;
  using AbstractValue::grad;

  /// \return lower bound of variable
  double max() const;
  /// \return upper bound of variable
  double min() const;

  /// \brief sets new upper bound for the variable
  /// \param new_max new upper bound for the variable
  void max(double new_max);
  /// \brief sets new lower bound for the variable
  /// \param new_min new lower bound for the variable
  void min(double new_min);
  /// \brief sets new bounds for the variable
  /// \param v1 first value, either minimum or maximum
  /// \param v2 second value, either minimum or maximum
  void bound(double v1, double v2);

  /// \brief sets current proxy variable so that nominal value equals new_val
  /// \param new_val new nominal value to set
  void val(double new_val) override;

  /// \returns transforms proxy variable into nominal value
  /// \param at_x proxy variable to be evaluated
  double val_at(double at_x) const override;

  /// \returns transforms proxy variable into gradient (slope)
  /// \param at_x proxy variable to be evaluated
  double grad_at(double at_x) const override;

  std::string to_string() const override;
  friend void to_json(nlohmann::json& j, const Value& s);
  friend void from_json(const nlohmann::json& j, Value& s);

 private:
  double max_{1.0};
  double min_{0.0};
};

class ValueGam : public AbstractValue
{
 public:
  ValueGam() = default;

  using AbstractValue::x;
  using AbstractValue::val;
  using AbstractValue::grad;

  void val(double new_val) override;
  double val_at(double at_x) const override;
  double grad_at(double at_x) const override;
};

class ValueBkg : public AbstractValue
{
 public:
  ValueBkg() = default;

  using AbstractValue::x;
  using AbstractValue::val;
  using AbstractValue::grad;

  void val(double new_val) override;
  double val_at(double at_x) const override;
  double grad_at(double at_x) const override;
};

}

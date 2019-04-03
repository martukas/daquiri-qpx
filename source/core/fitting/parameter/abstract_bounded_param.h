#pragma once

#include <core/fitting/parameter/abstract_param.h>

namespace DAQuiri
{

class AbstractBoundedParam : public AbstractParam
{
 public:
  AbstractBoundedParam() = default;

  using AbstractParam::x;
  using AbstractParam::val;
  using AbstractParam::grad;

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

  // \brief set bounds and value; makes comparisons and orders as needed
  void set(double v1, double v2, double v3);

  /// \returns if current value is at extermum
  /// \param min_epsilon maximum different between value and min to yield true
  /// \param max_epsilon maximum different between value and max to yield true
  bool at_extremum(double min_epsilon, double max_epsilon) const;

  std::string to_string() const override;
  friend void to_json(nlohmann::json& j, const AbstractBoundedParam& s);
  friend void from_json(const nlohmann::json& j, AbstractBoundedParam& s);

 private:
  double max_{1.0};
  double min_{0.0};
};

}

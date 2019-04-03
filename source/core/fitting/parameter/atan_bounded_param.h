#pragma once

#include <core/fitting/parameter/abstract_bounded_param.h>

namespace DAQuiri
{

class AtanBoundedParam : public AbstractBoundedParam
{
 public:
  double slope_ {1.0};

  AtanBoundedParam() = default;

  using AbstractParam::x;
  using AbstractParam::val;
  using AbstractParam::grad;
  using AbstractBoundedParam::min;
  using AbstractBoundedParam::max;
  using AbstractBoundedParam::bound;

  /// \brief sets current proxy variable so that nominal value equals new_val
  /// \param new_val new nominal value to set
  void val(double new_val) override;

  /// \returns transforms proxy variable into nominal value
  /// \param at_x proxy variable to be evaluated
  double val_at(double at_x) const override;

  /// \returns transforms proxy variable into gradient (slope)
  /// \param at_x proxy variable to be evaluated
  double grad_at(double at_x) const override;
};

}

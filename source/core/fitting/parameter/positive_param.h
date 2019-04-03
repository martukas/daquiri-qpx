#pragma once

#include <core/fitting/parameter/abstract_param.h>

namespace DAQuiri
{

class PositiveParam : public AbstractParam
{
 public:
  PositiveParam() = default;

  using AbstractParam::x;
  using AbstractParam::val;
  using AbstractParam::grad;

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

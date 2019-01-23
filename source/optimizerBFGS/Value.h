#pragma once

#include <cinttypes>
#include <vector>

namespace Hypermet
{

class AbstractValue
{
 public:
  AbstractValue() = default;
  virtual ~AbstractValue() = default;

  double uncert_value{0};
  int32_t x_index{-1};
  bool to_fit{true};

  double x() const;
  void x(double new_x);
  double val() const;
  double grad() const;

  virtual void val(double new_val) = 0;
  virtual double val_at(double at_x) const = 0;
  virtual double grad_at(double at_x) const = 0;

  void put(std::vector<double>& fit) const;
  void get(const std::vector<double>& fit);
  void get_uncert(const std::vector<double>& diagonals, double chisq_norm);

 private:
  double x_{0};
  double dx_{0};
};


class Value : public AbstractValue
{
 public:
  Value() = default;

  using AbstractValue::x;
  using AbstractValue::val;
  using AbstractValue::grad;

  double max() const;
  double min() const;

  void max(double new_max);
  void min(double new_min);
  void bound(double v1, double v2);

  void val(double new_val) override;
  double val_at(double at_x) const override;
  double grad_at(double at_x) const override;

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

struct PrecalcVals
{
  double width;
  double ampl;
  double half_ampl;
  double spread;
};

enum class Side {
  left,
  right
};

}
